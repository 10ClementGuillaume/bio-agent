"""Small offline checks for scoring, evidence extraction and citation provenance."""

import math
import io
import json
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from smolagents.utils import AgentError

import bio_agent as app
import bio_tools as sources
from bio_input import parse_task, scenario
from bio_console import DemoConsole
from rich.console import Console


class WorkflowTests(unittest.TestCase):
    def test_zero_hits_missing_data_and_experimental_channel(self):
        candidates = [{"gene": gene, "pubmed_query": gene, "experimental_score": score,
                       "combined_score": .999} for gene, score in [("ZERO", .742), ("FAILED", .99),
                                                                  ("TEXT_ONLY", 0), ("KNOWN", .994)]]
        searches = {"ZERO": {"status": "ok", "count": 0}, "FAILED": {"status": "error"},
                    "TEXT_ONLY": {"status": "ok", "count": 0}, "KNOWN": {"status": "ok", "count": 5}}
        rows = app.rank_candidates(candidates, searches)
        self.assertEqual([row["gene"] for row in rows], ["ZERO", "KNOWN", "TEXT_ONLY", "FAILED"])
        self.assertEqual(rows[0]["blindspot_score"], .742)
        self.assertAlmostEqual(rows[1]["blindspot_score"], .994 / (1 + math.log10(6)), places=5)
        self.assertFalse(rows[2]["eligible"])
        self.assertIsNone(rows[3]["pubmed_count"])
        self.assertEqual([row["gene"] for row in app.shortlist(rows)], ["ZERO", "KNOWN"])

    def test_pmc_table_boundaries_and_references(self):
        root = ET.fromstring('''<pmc-articleset><article><front><article-meta>
          <article-id pub-id-type="pmcid">PMC123</article-id><title-group><article-title>Study</article-title></title-group>
          </article-meta></front><body><p>SLC25A250 is a different token.</p>
          <table-wrap><table><tr><td>2.5</td><td>SLC25A25</td><td>expression</td></tr></table></table-wrap>
          </body><back><ref-list><p>SLC25A25 reference only</p></ref-list></back></article></pmc-articleset>''')
        result = sources.pmc_articles(root, "SLC25A25")[0]
        self.assertEqual(result["source_id"], "PMC123")
        self.assertEqual(result["snippets"], [{"kind": "table_row", "text": "2.5 | SLC25A25 | expression", "context": "", "truncated": False}])
        table = root.find(".//body/table-wrap")
        root.find(".//body").remove(table)
        ET.SubElement(root.find("article"), "floats-group").append(table)
        self.assertEqual(sources.pmc_articles(root, "SLC25A25")[0]["snippets"], result["snippets"])

    def test_offline_cache_miss_is_an_error(self):
        with tempfile.TemporaryDirectory() as directory:
            client = sources.EvidenceClient(directory, directory, offline=True)
            with self.assertRaisesRegex(OSError, "缓存缺失"):
                client.fetch("https://example.invalid/data")

    def test_agent_score_and_citation_checks(self):
        gene = "EXAMPLE"
        query = sources.pubmed_query([gene])
        candidate = {"gene": gene, "pubmed_query": query, "experimental_score": .8,
                     "source_id": "STRING:SLC7A11-EXAMPLE", "url": "https://string-db.org/"}
        calls = [
            {"tool": "get_interacting_proteins", "arguments": {"gene": "SLC7A11"},
             "result": {"status": "ok", "candidates": [candidate]}},
            {"tool": "get_pubmed_count_and_abstracts", "arguments": {"query": query},
             "result": {"status": "ok", "count": 0}},
            {"tool": "get_pubmed_count_and_abstracts", "arguments": {"query": sources.pubmed_query([gene], sources.BREAST)},
             "result": {"status": "ok", "count": 0}},
            {"tool": "get_gene_profile", "arguments": {"gene": gene}, "result": {"status": "error"}},
            {"tool": "search_pmc_evidence", "arguments": {"gene": gene}, "result": {"status": "ok", "articles": []}},
        ]
        with patch.object(sources, "client", SimpleNamespace(calls=calls, scenario=scenario("SLC7A11", "TNBC"))):
            agent = app.CodeAgent(tools=[], model=SimpleNamespace(model_id="test"),
                                  final_answer_checks=[app.check_screening])
            with self.assertRaisesRegex(AgentError, "分数有误"):
                agent._validate_final_answer([{"gene": gene, "pubmed_count": 0, "blindspot_score": .99}])
            agent._validate_final_answer([{"gene": gene, "pubmed_count": 0, "blindspot_score": .8}])
            agent.final_answer_checks = [app.check_review]
            with self.assertRaisesRegex(AgentError, "引用未出现"):
                agent._validate_final_answer({"reviews": [{"gene": gene, "assessment": "test", "gap": "test",
                    "evidence": [{"source_id": "PMID:invented", "kind": "mechanism", "claim": "test"}]}]})
            calls[1]["result"]["count"] = 5
            with self.assertRaisesRegex(AgentError, "低关注筛选条件"):
                agent._validate_final_answer({"reviews": [{"gene": gene, "assessment": "known", "gap": "test",
                    "evidence": []}], "hypothesis": {"gene": gene}})

    def test_scenario_reaches_string_candidate_queries(self):
        for gene, disease in [("GPX4", "胰腺癌"), ("EGFR", "Lung cancer")]:
            config = scenario(gene, disease)
            with tempfile.TemporaryDirectory() as directory:
                client = sources.EvidenceClient(directory, directory, scenario=config)
                responses = [[{"stable_address": "https://string-db.org"}], [
                    {"preferredName_B": "PARTNER", "stringId_B": "9606.test", "score": .9,
                     "escore": .8, "ascore": .1, "tscore": .4}]]
                with patch.object(sources, "client", client), patch.object(client, "fetch", side_effect=responses):
                    result = sources.get_interacting_proteins(gene)
                self.assertEqual(result["status"], "ok")
                self.assertEqual(result["seed"], gene)
                query = result["candidates"][0]["pubmed_query"]
                self.assertIn(config["disease_context"], query)
                self.assertNotIn("TNBC", query)
                self.assertEqual(client.calls[0]["arguments"], {"gene": gene})

    def test_natural_language_input_and_missing_context(self):
        model = SimpleNamespace(generate=lambda *args, **kwargs: SimpleNamespace(content=json.dumps(
            {"gene": "gpx4", "disease": "Pancreatic cancer", "broader_disease": None})))
        config = parse_task("围绕 GPX4 在胰腺癌中寻找候选", model)
        self.assertEqual(config["gene"], "GPX4")
        self.assertIn("Pancreatic Neoplasms", config["disease_context"])
        model.generate = lambda *args, **kwargs: SimpleNamespace(content='{"error":"请指定疾病"}')
        with self.assertRaisesRegex(ValueError, "请指定疾病"):
            parse_task("研究 GPX4", model)
        with self.assertRaises(ValueError):
            scenario("GPX4", 'cancer" OR all[All Fields]')

    def test_readable_terminal_ranking(self):
        stream = io.StringIO()
        ui = DemoConsole(Console(file=stream, width=100, color_system=None))
        rows = app.rank_candidates([{"gene": "EXAMPLE", "pubmed_query": "query", "experimental_score": .8}],
                                   {"query": {"status": "ok", "count": 0}})
        ui.header()
        ui.task(scenario("GPX4", "胰腺癌"), False)
        ui.ranking(rows)
        rendered = stream.getvalue()
        for expected in ("GPX4", "Pancreatic cancer", "低关注候选", "0.800", "PubMed 篇数"):
            self.assertIn(expected, rendered)
        self.assertNotIn("experimental_score", rendered)
        ui.finish(rows, {"reviews": [], "hypothesis": {
            "gene": None, "statement": "本次无低关注候选", "prediction": "旧预测不应显示",
            "validation": "旧实验不应显示", "limitations": "等待更多证据", "source_ids": []}},
            Path.cwd() / "outputs" / "test", [])
        self.assertIn("本次无低关注候选", stream.getvalue())
        self.assertNotIn("旧预测不应显示", stream.getvalue())


if __name__ == "__main__":
    unittest.main()
