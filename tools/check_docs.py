"""Check local Markdown links, documented module paths and Python syntax."""

from __future__ import annotations

import ast
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
FENCE = re.compile(r"^```([^\n]*)\n(.*?)^```\s*$", re.MULTILINE | re.DOTALL)


def python_blocks(text):
    return [
        body for language, body in FENCE.findall(text) if language.strip() == "python"
    ]


def anchors(text):
    seen, result = {}, set()
    for title in re.findall(
        r"^#{1,6}\s+(.+?)\s*#*\s*$", FENCE.sub("", text), re.MULTILINE
    ):
        slug = re.sub(r"[^\w\- ]", "", title.lower()).replace(" ", "-")
        count = seen.get(slug, 0)
        seen[slug] = count + 1
        result.add(f"{slug}-{count}" if count else slug)
    result.update(re.findall(r"\bid=[\"\']([^\"\']+)", text))
    return result


def check(root=ROOT):
    files = [
        root / name
        for name in (
            "README.md",
            "CONTEXT.md",
            "CLAUDE.md",
            "tools/README.md",
            "tests/AGENTS.md",
            "tests/browser/README.md",
        )
    ]
    files.extend(sorted((root / "docs").rglob("*.md")))
    errors, links, snippets = [], 0, 0
    for path in files:
        text = path.read_text(encoding="utf-8")
        for code in python_blocks(text):
            snippets += 1
            try:
                ast.parse(
                    code, filename=str(path.relative_to(root)), feature_version=(3, 9)
                )
            except SyntaxError as error:
                errors.append(f"{path.relative_to(root)}: {error}")
        prose = re.sub(r"`[^`\n]+`", "", FENCE.sub("", text))
        for target in re.findall(r"\]\(([^\s)]+)(?:\s+[^)]*)?\)", prose):
            link = urlsplit(target.strip("<>"))
            if link.scheme or link.netloc:
                continue
            links += 1
            destination = (
                (path.parent / unquote(link.path)).resolve() if link.path else path
            )
            if not destination.exists():
                errors.append(f"{path.relative_to(root)}: missing {target}")
            elif link.fragment and destination.suffix == ".md":
                if unquote(link.fragment) not in anchors(
                    destination.read_text(encoding="utf-8")
                ):
                    errors.append(f"{path.relative_to(root)}: missing anchor {target}")
    architecture = (root / "docs/architecture.md").read_text(encoding="utf-8")
    for name in re.findall(
        r"\| `((?:pyPept\.[\w.]+)|(?:[\w-]+\.js))` \|", architecture
    ):
        path = root / (
            "src/" + name.replace(".", "/")
            if name.startswith("pyPept.")
            else "src/pyPept/web/static/" + name
        )
        if not (path.exists() or path.with_suffix(".py").exists()):
            errors.append(f"docs/architecture.md: missing module {name}")
    if errors:
        raise SystemExit("\n".join(errors))
    print(
        f"Checked {len(files)} documents, {links} local links and {snippets} Python snippets."
    )


if __name__ == "__main__":
    check()
