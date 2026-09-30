"""Save evidence only for marked fuzz tests; wheel tests need no Hypothesis."""

import pytest


def pytest_runtest_setup(item):
    if item.get_closest_marker("fuzz"):
        from _fuzzing import start

        start(item.nodeid, item.config.invocation_params.args)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    report = (yield).get_result()
    if report.when == "call" and item.get_closest_marker("fuzz"):
        from _fuzzing import finish

        finish(item.nodeid, report, item.config.invocation_params.args)
