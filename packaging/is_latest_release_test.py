#!/usr/bin/env python3
"""Local tests for packaging/is_latest_release.py.

    python3 packaging/is_latest_release_test.py

Covers the `latest` decision of the release workflow's Docker job:
newest stable -> latest, older stable (backport) -> no latest, pre-release ->
no latest, first C++ release among the Julia-era v0.* releases -> latest,
plus drafts/prereleases being ignored and the command-line interface.
"""

import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "is_latest_release", HERE / "is_latest_release.py")
is_latest_release = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(is_latest_release)


def rel(tag, prerelease=False, draft=False):
    """One entry as `gh release list --json tagName,isDraft,isPrerelease`."""
    return {"tagName": tag, "isPrerelease": prerelease, "isDraft": draft}


class VersionKeyTest(unittest.TestCase):
    def test_plain_versions(self):
        self.assertEqual(is_latest_release.version_key("0.8.3"),
                         ((0, 8, 3), True, ""))
        # different lengths compare as versions, not as strings
        self.assertLess(is_latest_release.version_key("0.4"),
                        is_latest_release.version_key("0.8.3"))
        self.assertLess(is_latest_release.version_key("0.8.2"),
                        is_latest_release.version_key("0.8.10"))

    def test_suffix_sorts_before_plain(self):
        self.assertLess(is_latest_release.version_key("0.9.0-rc1"),
                        is_latest_release.version_key("0.9.0"))

    def test_unparsable(self):
        self.assertIsNone(is_latest_release.version_key("nightly"))

    def test_tag_forms(self):
        for tag in ("cpp-0.8.3", "cpp-v0.8.3", "v0.8.3", "0.8.3"):
            self.assertEqual(is_latest_release.version_of_tag(tag),
                             ((0, 8, 3), True, ""), tag)


class IsLatestTest(unittest.TestCase):
    def test_newest_stable_is_latest(self):
        releases = [rel("cpp-0.9.1", prerelease=True),  # pre-release: ignored
                    rel("cpp-0.9.0"), rel("cpp-0.8.3"), rel("v0.7.1")]
        self.assertTrue(is_latest_release.is_latest("0.9.0", False, releases))

    def test_older_stable_backport_is_not_latest(self):
        releases = [rel("cpp-0.9.0"), rel("cpp-0.8.3"), rel("cpp-0.8.2")]
        self.assertFalse(is_latest_release.is_latest("0.8.1", False, releases))

    def test_prerelease_is_not_latest(self):
        # newest of all, but a GitHub pre-release
        releases = [rel("cpp-0.9.0"), rel("v0.8.3", prerelease=True)]
        self.assertFalse(is_latest_release.is_latest("0.9.1", True, releases))
        # the pre-release itself, present in its own release list
        self.assertFalse(is_latest_release.is_latest("0.8.3", True, releases))

    def test_first_cpp_release_among_julia_releases_is_latest(self):
        releases = [rel("v0.7.1"), rel("v0.7.0"), rel("v0.6.2"),
                    rel("v0.5.0"), rel("v0.4.0"), rel("v0.2.0")]
        self.assertTrue(is_latest_release.is_latest("0.8.0", False, releases))
        # a Julia-era version would not be latest among the Julia releases
        self.assertFalse(is_latest_release.is_latest("0.7.0", False, releases))

    def test_drafts_are_ignored(self):
        releases = [rel("cpp-1.0.0", draft=True),
                    rel("cpp-1.1.0", prerelease=True)]
        self.assertTrue(is_latest_release.is_latest("0.9.0", False, releases))

    def test_candidate_equal_to_newest_stable_is_latest(self):
        # the workflow lists releases including the one being published
        releases = [rel("cpp-0.8.3"), rel("v0.7.1")]
        self.assertTrue(is_latest_release.is_latest("0.8.3", False, releases))

    def test_unparsable_version_raises(self):
        with self.assertRaises(ValueError):
            is_latest_release.is_latest("nightly", False, [])


class CommandLineTest(unittest.TestCase):
    def run_main(self, releases, *argv):
        with tempfile.NamedTemporaryFile("w", suffix=".json",
                                         delete=False) as f:
            json.dump(releases, f)
            path = f.name
        out = subprocess.run(
            [sys.executable, str(HERE / "is_latest_release.py"),
             *argv, "--releases-file", path],
            capture_output=True, text=True)
        pathlib.Path(path).unlink()
        self.assertEqual(out.returncode, 0, out.stderr)
        return out.stdout.strip()

    def test_prints_true_and_false(self):
        releases = [rel("cpp-0.8.3"), rel("v0.7.1")]
        self.assertEqual(self.run_main(releases, "0.8.4"), "true")
        self.assertEqual(self.run_main(releases, "0.8.2"), "false")
        self.assertEqual(self.run_main(releases, "0.8.4",
                                       "--prerelease"), "false")


if __name__ == "__main__":
    unittest.main(verbosity=2)
