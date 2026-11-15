"""Compact terminal presentation for a recorded demo."""

import re
import time
from pathlib import Path

from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text


def excerpt(value, limit):
    return value if len(value) <= limit else value[:limit].rstrip() + "…"


class DemoConsole:
    def __init__(self, console=None):
        self.console = console or Console(highlight=False)
        self.started = time.monotonic()
        self.status = self.console.status("准备任务…", spinner="dots")
        self.total = 0
        self.counted = 0
        self.genes = []
        self.stage = 0

    def header(self):
        self.console.print(Panel("关联蛋白筛选 → 文献核查 → 待验证假说", title="BioBlindSpot Lite",
                                 border_style="cyan", padding=(1, 2)))

    def task(self, config, data_only):
        text = Text()
        if config["task"]:
            text.append(f"研究问题  {config['task']}\n")
        text.append(f"核心蛋白  {config['gene']}\n疾病场景  {config['disease']}\n")
        text.append(f"运行模式  {'数据基线' if data_only else 'CodeAgent + DeepSeek'} · 人类蛋白\n")
        text.append("筛选条件  实验评分 ≥ 0.7，初筛 PubMed 命中 ≤ 2")
        self.console.print(Panel(text, title="任务输入", border_style="cyan"))

    def phase(self, number, label):
        self.stage = number
        self.console.rule(f"{number}/3  {label}", style="cyan")
        self.status.update(f"{label} · 正在处理…")

    def tool_event(self, name, arguments, result=None):
        label = arguments.get("gene")
        if not label:
            query = arguments.get("query", "")
            label = next((gene for gene in self.genes if re.search(
                rf"(?<![A-Za-z0-9]){re.escape(gene)}(?![A-Za-z0-9])", query, re.I)), None)
            if not label:
                match = re.search(r'"([^\"]+)"', query)
                label = match.group(1) if match else "文献"
        source = {"get_interacting_proteins": "STRING", "get_pubmed_count_and_abstracts": "PubMed",
                  "get_gene_profile": "UniProt", "search_pmc_evidence": "PMC"}[name]
        if result is None:
            self.status.update(f"{source} · {label} · 正在检索…")
            return
        if result["status"] != "ok":
            self.console.print(Text(f"  ! {source} · {label}：检索失败，详情已保存到工具日志。", style="yellow"))
        elif name == "get_interacting_proteins":
            self.total = len(result["candidates"])
            self.genes = [row["gene"] for row in result["candidates"]]
            self.console.print(f"  STRING 返回 {self.total} 个关联蛋白，开始逐一统计文献。")
        elif self.stage == 1 and name == "get_pubmed_count_and_abstracts":
            self.counted += 1
        elif self.stage == 2:
            detail = (f"命中 {result['count']} 篇，读取 {len(result.get('articles', []))} 篇"
                      if "count" in result else "功能与别名已获取" if result.get("found") else "未找到功能注释")
            self.console.print(Text(f"  {source} · {label}：{detail}"))
        progress = f"文献统计 {self.counted}/{self.total} · " if self.stage == 1 and self.total else ""
        self.status.update(f"{progress}Agent 正在分析证据…")

    def ranking(self, rows):
        table = Table(title="候选排名 · 按 BlindSpot Score 降序", box=box.ROUNDED, header_style="bold cyan")
        for label in ("蛋白", "实验评分", "PubMed 篇数", "盲区分数", "筛选结果"):
            table.add_column(label, justify="left" if label in {"蛋白", "筛选结果"} else "right")
        for row in rows:
            label = "低关注候选" if row["low_attention"] else "已有文献" if row["eligible"] else "实验评分不足"
            if row["pubmed_count"] is None:
                label = "检索失败"
            table.add_row(row["gene"], f"{row['experimental_score']:.3f}",
                          str(row["pubmed_count"]) if row["pubmed_count"] is not None else "—",
                          f"{row['blindspot_score']:.3f}" if row["blindspot_score"] is not None else "—",
                          label, style="yellow" if row["low_attention"] else None)
        self.console.print(table)
        self.console.print("  盲区分数 = 实验评分 ÷ (1 + log10(PubMed 篇数 + 1))", style="dim")
        self.console.print("  低关注候选仍需核查摘要和全文中的已有研究。", style="dim")

    def finish(self, rows, review, output, calls):
        if review:
            self.console.print("终端展示证据摘要，完整论述及来源链接见 report.md。", style="dim")
            for item in review["reviews"]:
                text = Text(excerpt(item["assessment"], 180) + "\n\n待核查：" + excerpt(item["gap"], 100))
                for evidence in item["evidence"][:2]:
                    text.append(f"\n• {evidence['source_id']}：{excerpt(evidence['claim'], 90)}", style="dim")
                self.console.print(Panel(text, title=f"证据核查 · {item['gene']}", border_style="blue"))
            card = review["hypothesis"]
            text = Text()
            for key, label in (("statement", "假说"), ("prediction", "预测"),
                               ("validation", "验证"), ("limitations", "证据缺口")):
                if card["gene"] is None and key in {"prediction", "validation"}:
                    continue
                if card["gene"] is None and key == "statement":
                    label = "筛选结论"
                text.append(label + "  ", style="bold cyan")
                text.append(excerpt(str(card[key]), 260) + "\n\n")
            text.append("依据  " + "、".join(card["source_ids"]), style="dim")
            self.console.print(Panel(text, title=f"待验证假说 · {card['gene'] or '暂无合适候选'}", border_style="cyan"))
        failed = sum(call["result"]["status"] != "ok" for call in calls)
        summary = (f"完成 · {len(rows)} 个关联蛋白 · {sum(r['low_attention'] for r in rows)} 个低关注候选"
                   f" · 用时 {time.monotonic() - self.started:.0f} 秒")
        if failed:
            summary += f"\n{failed} 次工具请求失败，详见 tool_calls.jsonl。"
        directory = output.relative_to(Path.cwd()) if output.is_relative_to(Path.cwd()) else output
        summary += (f"\n\n输出目录  {directory}\n报告  report.md    散点图  blindspots.png"
                    "\n候选表  candidates.csv    完整证据  result.json\n执行日志  run.log    工具记录  tool_calls.jsonl")
        self.console.print(Panel(Text(summary), title="运行结果", border_style="green"))
