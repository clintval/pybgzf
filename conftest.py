from sybil import Sybil
from sybil.parsers.markdown import PythonCodeBlockParser

pytest_collect_file = Sybil(parsers=[PythonCodeBlockParser()], patterns=["*.md"]).pytest()
