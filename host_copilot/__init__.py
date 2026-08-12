"""Host Copilot public Python API."""

__version__ = "0.1.0"

from .catalog import GalaxyFinder
from .image import Imager
from .models import (
    CircleLocalization,
    EllipseLocalization,
    HostCandidate,
    HostSearchResult,
    ProviderStatus,
    SearchConfig,
    TransientContext,
)
from .pipeline import HostPipeline

__all__ = [
    "CircleLocalization",
    "EllipseLocalization",
    "GalaxyFinder",
    "HostCandidate",
    "HostPipeline",
    "HostSearchResult",
    "Imager",
    "ProviderStatus",
    "SearchConfig",
    "TransientContext",
]
