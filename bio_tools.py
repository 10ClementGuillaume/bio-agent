"""Public evidence sources for a human gene and disease scenario."""

import hashlib
import inspect
import json
import re
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from functools import wraps
from pathlib import Path
from urllib.parse import urlencode, urlparse

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from smolagents import tool

TNBC = '("Triple Negative Breast Neoplasms"[MeSH Terms] OR "triple negative breast cancer"[Title/Abstract] OR TNBC[Title/Abstract])'
BREAST = '("Breast Neoplasms"[MeSH Terms] OR "breast cancer"[Title/Abstract] OR "breast carcinoma"[Title/Abstract])'
PMC_CONTEXT = '("triple negative breast cancer" OR TNBC OR "MDA-MB-231")'
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


class EvidenceClient:
    """One serial HTTP session with raw-response caching and provenance."""

    def __init__(self, cache_dir, output_dir, offline=False, scenario=None, on_tool=None):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.log_path = Path(output_dir) / "tool_calls.jsonl"
        self.offline = offline
        self.scenario = scenario or {"disease_context": TNBC, "broad_context": BREAST}
        self.on_tool = on_tool
        self.calls = []
        self.requests = []
        self.last_request = {}
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "BioBlindSpot-Lite/0.1 (course prototype)"
        retry = Retry(total=2, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
        self.session.mount("https://", HTTPAdapter(max_retries=retry))

    def fetch(self, url, params=None, xml=False):
        params = dict(params or {})
        if url.startswith(EUTILS):
            params["tool"] = "BioBlindSpotLite"
        full_url = url + ("?" + urlencode(sorted(params.items())) if params else "")
        cache_key = hashlib.sha256(full_url.encode()).hexdigest()
        cache_path = self.cache_dir / f"{cache_key}.json"
        cached = cache_path.exists()
        if cached:
            record = json.loads(cache_path.read_text())
        else:
            if self.offline:
                raise OSError(f"缓存缺失：{full_url}")
            host = urlparse(url).netloc
            interval = 0.4 if host == "eutils.ncbi.nlm.nih.gov" else 1.0
            time.sleep(max(0, interval - (time.monotonic() - self.last_request.get(host, 0))))
            self.last_request[host] = time.monotonic()
            response = self.session.get(url, params=params, timeout=(10, 45))
            response.raise_for_status()
            record = {"url": full_url, "retrieved_at": datetime.now(timezone.utc).isoformat(),
                      "body": response.text}
        parsed = ET.fromstring(record["body"]) if xml else json.loads(record["body"])
        if xml:
            error = parsed.find(".//ERROR")
            if error is not None:
                raise ValueError(error.text)
        elif isinstance(parsed, dict):
            error = parsed.get("error") or parsed.get("esearchresult", {}).get("errorlist")
            if error:
                raise ValueError(str(error))
        if not cached:
            cache_path.write_text(json.dumps(record, ensure_ascii=False), encoding="utf-8")
        self.requests.append({k: record[k] for k in ("url", "retrieved_at")} |
                             {"cache_key": cache_key, "cached": cached})
        return parsed


client = None


def recorded(function):
    @wraps(function)
    def call(*args, **kwargs):
        bound = inspect.signature(function).bind(*args, **kwargs)
        bound.apply_defaults()
        start = len(client.requests)
        if client.on_tool:
            client.on_tool(function.__name__, bound.arguments)
        try:
            result = {"status": "ok", **function(*args, **kwargs)}
        except (requests.RequestException, ValueError, KeyError, ET.ParseError, OSError) as error:
            result = {"status": "error", "error": str(error)}
        entry = {"tool": function.__name__, "arguments": bound.arguments,
                 "result": result, "sources": client.requests[start:]}
        client.calls.append(entry)
        with client.log_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(entry, ensure_ascii=False) + "\n")
        if client.on_tool:
            client.on_tool(function.__name__, bound.arguments, result)
        return result
    return call


def text_of(element):
    return " ".join("".join(element.itertext()).split()) if element is not None else ""


def pubmed_query(names, context=TNBC):
    terms = " OR ".join(f'"{name}"[Title/Abstract]' for name in names)
    return f"({terms}) AND {context}"


@tool
@recorded
def get_interacting_proteins(gene: str) -> dict:
    """Get up to 20 human STRING functional neighbors with combined score >= 0.7.

    Args:
        gene: Human seed gene symbol, for example SLC7A11 or GPX4.
    """
    version = client.fetch("https://string-db.org/api/json/version")[0]
    base = version["stable_address"].rstrip("/")
    data = client.fetch(base + "/api/json/interaction_partners", {
        "identifiers": gene, "species": 9606, "required_score": 700, "limit": 20,
        "caller_identity": "BioBlindSpotLite"})
    rows = [{"gene": row["preferredName_B"], "string_id": row["stringId_B"],
             "source_id": f"STRING:{gene}-{row['preferredName_B']}",
             "url": base + "/api/json/network?" + urlencode({"identifiers": gene + "\r" + row["preferredName_B"], "species": 9606}),
             "pubmed_query": pubmed_query([row["preferredName_B"]], client.scenario["disease_context"]),
             "combined_score": row["score"], "experimental_score": row["escore"],
             "coexpression_score": row["ascore"], "textmining_score": row["tscore"]}
            for row in data]
    return {"seed": gene, "species": 9606, "string_version": version, "candidates": rows}


@tool
@recorded
def get_pubmed_count_and_abstracts(query: str, max_abstracts: int = 0) -> dict:
    """Count PubMed matches; optionally read up to three abstracts. Failures have no count.

    Args:
        query: Complete PubMed query including field tags.
        max_abstracts: Number of abstracts to fetch, from zero to three.
    """
    limit = min(max(max_abstracts, 0), 3)
    search = client.fetch(EUTILS + "esearch.fcgi", {
        "db": "pubmed", "term": query, "retmode": "json", "retmax": limit,
        "sort": "relevance"})["esearchresult"]
    articles = []
    if search["idlist"]:
        root = client.fetch(EUTILS + "efetch.fcgi", {
            "db": "pubmed", "id": ",".join(search["idlist"]), "retmode": "xml"}, xml=True)
        for article in root.findall(".//PubmedArticle"):
            pmid = text_of(article.find(".//MedlineCitation/PMID"))
            articles.append({"source_id": "PMID:" + pmid, "pmid": pmid,
                             "title": text_of(article.find(".//ArticleTitle")),
                             "abstract": "\n".join(text_of(p) for p in article.findall(".//AbstractText")),
                             "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
    return {"query": query, "query_translation": search["querytranslation"],
            "warnings": search.get("warninglist", {}), "count": int(search["count"]),
            "articles": articles}


@tool
@recorded
def get_gene_profile(gene: str) -> dict:
    """Read a reviewed human UniProt entry: names, function, location and evidence codes.

    Args:
        gene: Exact human gene symbol.
    """
    data = client.fetch("https://rest.uniprot.org/uniprotkb/search", {
        "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true", "format": "json",
        "fields": "accession,gene_names,protein_name,cc_function,cc_subcellular_location", "size": 1})
    if not data["results"]:
        return {"gene": gene, "found": False, "aliases": []}
    entry = data["results"][0]
    accession = entry["primaryAccession"]
    aliases = [item["value"] for names in entry.get("genes", []) for item in names.get("synonyms", [])]
    return {"gene": gene, "found": True, "accession": accession, "source_id": "UniProt:" + accession,
            "url": f"https://www.uniprot.org/uniprotkb/{accession}/entry", "aliases": aliases,
            "protein": entry.get("proteinDescription", {}), "annotations": entry.get("comments", [])}


def pmc_articles(root, gene):
    """Return gene-containing body paragraphs and table rows, excluding references."""
    pattern = re.compile(r"(?<![A-Za-z0-9])" + re.escape(gene) + r"(?![A-Za-z0-9])", re.I)
    articles = []
    for article in root.findall(".//article"):
        parents = {child: parent for parent in article.iter() for child in parent}
        ids = {node.get("pub-id-type"): text_of(node) for node in article.findall("./front/article-meta/article-id")}
        pmcid = ids.get("pmcid", ids.get("pmc", ""))
        if pmcid and not pmcid.startswith("PMC"):
            pmcid = "PMC" + pmcid
        snippets = []
        for path, kind in [(".//body//tr", "table_row"), ("./floats-group//tr", "table_row"),
                           (".//body//p", "paragraph")]:
            for node in article.findall(path):
                if kind == "paragraph" and node.find(".//table-wrap") is not None:
                    continue
                # Spaces between cells keep gene symbols separate from numeric values.
                text = " | ".join(text_of(cell) for cell in node) if kind == "table_row" else text_of(node)
                match = pattern.search(text)
                if match:
                    ancestor = node
                    while ancestor in parents and ancestor.tag not in ("sec", "table-wrap"):
                        ancestor = parents[ancestor]
                    context = text_of(ancestor.find("caption")) if ancestor.tag == "table-wrap" else text_of(ancestor.find("title"))
                    left, right = max(0, match.start() - 250), match.end() + 650
                    snippets.append({"kind": kind, "text": text[left:right],
                                     "context": context,
                                     "truncated": left > 0 or right < len(text)})
        snippets.sort(key=lambda item: 0 if re.search(r"results|discussion", item["context"], re.I) else 1)
        articles.append({"source_id": pmcid, "pmid": ids.get("pmid"),
                         "title": text_of(article.find("./front/article-meta/title-group/article-title")),
                         "url": f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/", "snippets": snippets[:6]})
    return articles


@tool
@recorded
def search_pmc_evidence(query: str, gene: str) -> dict:
    """Search PMC and extract gene mentions from at most three matching full texts.

    Args:
        query: Full-text search query, usually gene AND the selected disease context.
        gene: Gene symbol whose body paragraphs and table rows should be extracted.
    """
    search = client.fetch(EUTILS + "esearch.fcgi", {
        "db": "pmc", "term": query, "retmode": "json", "retmax": 3, "sort": "relevance"})["esearchresult"]
    articles = []
    if search["idlist"]:
        root = client.fetch(EUTILS + "efetch.fcgi", {
            "db": "pmc", "id": ",".join(search["idlist"]), "retmode": "xml"}, xml=True)
        articles = pmc_articles(root, gene)
    return {"gene": gene, "query": query, "query_translation": search["querytranslation"],
            "warnings": search.get("warninglist", {}), "count": int(search["count"]),
            "articles": articles, "read_limit": 3}


TOOLS = [get_interacting_proteins, get_pubmed_count_and_abstracts, get_gene_profile, search_pmc_evidence]
