import os

from hypothesis import Phase
from hypothesis import settings

settings.register_profile(
    "ci", derandomize=True, phases=[Phase.explicit, Phase.reuse, Phase.generate]
)
settings.load_profile("ci" if os.environ.get("CI") else "default")
