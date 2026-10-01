#!/usr/bin/env python3
"""Tests for packaging/is_latest_release.py.

    python3 packaging/is_latest_release_test.py
"""

import contextlib
import importlib.util
import io
import pathlib
import unittest
from unittest import mock

_spec = importlib.util.spec_from_file_location(
    "is_latest_release", pathlib.Path(__file__).resolve().parent / "is_latest_release.py")
ilr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ilr)


def rel(tag, prerelease=False, draft=False):
    """One entry of `gh release list --json tagName,isDraft,isPrerelease`."""
    return {"tagName": tag, "isPrerelease": prerelease, "isDraft": draft}


JULIA = [rel(t) for t in ("v0.7.1", "v0.7.0", "v0.6.2", "v0.5.0", "v0.4.0", "v0.2.0")]


class VersionKeyTest(unittest.TestCase):
    def test_parsing(self):
        self.assertEqual(ilr.version_key("0.8.3"), ((0, 8, 3), True, ""))
        self.assertIsNone(ilr.version_key("nightly"))
        for tag in ("cpp-0.8.3", "cpp-v0.8.3", "v0.8.3", "0.8.3"):
            self.assertEqual(ilr.version_of_tag(tag), ((0, 8, 3), True, ""), tag)

    def test_order(self):
        # numeric segments of any length; a suffix sorts before the plain version
        for older, newer in (("0.4", "0.8.3"), ("0.8.2", "0.8.10"), ("0.9.0-rc1", "0.9.0")):
            self.assertLess(ilr.version_key(older), ilr.version_key(newer))


class IsLatestTest(unittest.TestCase):
    def test_rule(self):
        cases = [  # (what, version, prerelease, releases, expected)
            ("newest stable, pre-releases ignored", "0.9.0", False,
             [rel("cpp-0.9.1", prerelease=True), rel("cpp-0.9.0"), rel("cpp-0.8.3"), rel("v0.7.1")], True),
            ("backport", "0.8.1", False, [rel("cpp-0.9.0"), rel("cpp-0.8.3"), rel("cpp-0.8.2")], False),
            ("newest, but a pre-release", "0.9.1", True, [rel("cpp-0.9.0"), rel("v0.8.3", prerelease=True)], False),
            ("pre-release in its own list", "0.8.3", True, [rel("cpp-0.9.0"), rel("v0.8.3", prerelease=True)], False),
            ("first C++ release after the Julia ones", "0.8.0", False, JULIA, True),
            ("Julia-era version", "0.7.0", False, JULIA, False),
            ("drafts ignored", "0.9.0", False, [rel("cpp-1.0.0", draft=True), rel("cpp-1.1.0", prerelease=True)], True),
            ("equal to the newest stable", "0.8.3", False, [rel("cpp-0.8.3"), rel("v0.7.1")], True),
        ]
        for what, version, prerelease, releases, expected in cases:
            with self.subTest(what):
                self.assertIs(ilr.is_latest(version, prerelease, releases), expected)

    def test_unparsable_version_raises(self):
        with self.assertRaises(ValueError):
            ilr.is_latest("nightly", False, [])


class CommandLineTest(unittest.TestCase):
    def test_prints_true_and_false(self):
        releases = [rel("cpp-0.8.3"), rel("v0.7.1")]
        for argv, expected in ((["0.8.4"], "true"), (["0.8.2"], "false"), (["0.8.4", "--prerelease"], "false")):
            with mock.patch.object(ilr, "fetch_releases", return_value=releases) as fetch, \
                    contextlib.redirect_stdout(io.StringIO()) as out:
                ilr.main([*argv, "--repo", "owner/repo"])
            self.assertEqual(out.getvalue(), expected + "\n")
            fetch.assert_called_once_with("owner/repo")


if __name__ == "__main__":
    unittest.main(verbosity=2)
