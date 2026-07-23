from pathlib import Path
from typing import List as _List
from typing import Mapping, Optional

__all__: _List[str]

def path(resource: Optional[str] = ...) -> Path: ...
def clear(*resources: str) -> None: ...
def __dir__() -> _List[str]: ...
def _cached_download(
    url: str,
    *,
    resource: str,
    category: str,
    max_age: Optional[float],
    suffix: Optional[str] = ...,
) -> Path: ...
def _platform_cache_path(
    *,
    platform: str,
    environment: Mapping[str, str],
    home: Path,
) -> Path: ...
def _discard_cached_download(cache_file: Path) -> None: ...
def _unlink(file: Path) -> None: ...
