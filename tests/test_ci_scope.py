"""CI skips expensive jobs only with a complete, documentation-only diff."""

import subprocess

import pytest

from tools.ci_scope import changed_paths, classify


@pytest.mark.parametrize(
    "paths,application,documentation",
    [
        (["README.md", "docs/nested/guide.md"], False, True),
        (["docs/guide.md", "src/pyPept/editor.py"], True, True),
        (["docs/example.py"], True, False),
        (["src/pyPept/data/README.md"], True, False),
        ([".github/workflows/tests.yml"], True, False),
        ([], False, False),
        (None, True, True),
    ],
)
def test_only_known_documentation_paths_skip_application_jobs(
    paths, application, documentation
):
    assert classify(paths) == {
        "application": application,
        "documentation": documentation,
    }


def test_git_comparison_includes_deleted_code_when_renamed_into_docs(tmp_path):
    def git(*args):
        return subprocess.check_output(
            ["git", "-C", str(tmp_path), *args], text=True
        ).strip()

    git("init", "-q")
    git("config", "user.name", "CI test")
    git("config", "user.email", "ci@example.invalid")
    (tmp_path / "source.py").write_text('print("example")\n')
    git("add", ".")
    git("commit", "-qm", "base")
    before = git("rev-parse", "HEAD")
    git("mv", "source.py", "README.md")
    git("commit", "-qm", "move")
    after = git("rev-parse", "HEAD")
    for event, data in [
        ("push", {"before": before, "after": after}),
        (
            "pull_request",
            {"pull_request": {"base": {"sha": before}, "head": {"sha": after}}},
        ),
    ]:
        paths = changed_paths(event, data, cwd=tmp_path)
        assert set(paths) == {"source.py", "README.md"}
        assert classify(paths) == {"application": True, "documentation": True}
    for event, data in [
        ("push", {"before": "0" * 40, "after": after}),
        ("push", {"before": "1" * 40, "after": after}),
        ("workflow_dispatch", {}),
        ("schedule", {}),
    ]:
        assert classify(changed_paths(event, data, cwd=tmp_path)) == {
            "application": True,
            "documentation": True,
        }
