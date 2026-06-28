#!/usr/bin/env python3
"""Submit a manual NIST kinetics reaction search and save the raw HTML."""

from __future__ import annotations

import argparse
import html
import http.cookiejar
import re
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE_URL = "https://kinetics.nist.gov/kinetics/"
USER_AGENT = "VULCAN-Reaction-Search/1.0 manual-research"


def default_temp_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "data_temp"


def slug(parts: list[str]) -> str:
    value = "_".join(parts)
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", value).strip("_") or "reaction"


def strip_tags(value: str) -> str:
    value = re.sub(r"<[^>]+>", "", value)
    value = html.unescape(value).replace("\u2192", "->").replace("\u00b7", "")
    return " ".join(value.split())


def build_form(reactants: list[str], products: list[str]) -> dict[str, str]:
    fields = [("reactants", item) for item in reactants]
    fields += [("products", item) for item in products]
    if not 1 <= len(fields) <= 4:
        raise ValueError("The NIST quick-search example supports 1-4 total species")

    form = {"database": "kinetics", "numberOfFields": str(len(fields))}
    for index, (field, value) in enumerate(fields, start=1):
        form[f"boolean{index}"] = " " if index == 1 else "and"
        form[f"lp{index}"] = " "
        form[f"field{index}"] = field
        form[f"relate{index}"] = "="
        form[f"rp{index}"] = " "
        form[f"text{index}"] = value
    return form


def matched_groups(page: str) -> list[tuple[str, str]]:
    pattern = re.compile(
        r'<a href="(?P<url>[^"]*ReactionSearch[^"]*)">(?P<count>[^<]*records? matched)</a>'
        r"</td>\s*<td[^>]*>.*?</td>\s*<td[^>]*>(?P<reaction>.*?)</td>",
        re.IGNORECASE | re.DOTALL,
    )
    return [
        (strip_tags(match.group("count")), strip_tags(match.group("reaction")))
        for match in pattern.finditer(page)
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reactant", action="append", required=True)
    parser.add_argument("--product", action="append", default=[])
    parser.add_argument("--output-dir", type=Path, default=default_temp_dir())
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cookie_jar = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cookie_jar))
    request = urllib.request.Request(
        urllib.parse.urljoin(BASE_URL, "Search.jsp"),
        data=urllib.parse.urlencode(build_form(args.reactant, args.product)).encode(),
        headers={"User-Agent": USER_AGENT},
    )
    with opener.open(request, timeout=60) as response:
        page = response.read().decode("utf-8", errors="replace")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    name = slug(args.reactant + ["to"] + args.product)
    output = args.output_dir / f"nist_{name}_{stamp}.html"
    output.write_text(page, encoding="utf-8")

    print(f"Saved raw NIST result: {output}")
    groups = matched_groups(page)
    if not groups:
        print("No reaction groups parsed. Inspect the saved HTML for site changes or zero results.")
        return
    for count, reaction in groups:
        print(f"{count}: {reaction}")
    print("Open a matching group in NIST and inspect individual records before selection.")


if __name__ == "__main__":
    main()
