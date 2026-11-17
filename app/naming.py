"""Rules for user-supplied names that are turned into filesystem paths.

A lane name becomes exactly one directory directly under the published-lane
root.  Only an ordinary single path component is acceptable: ``.``, ``..``
and anything containing a separator would resolve to the lanes root itself
or outside it, so they are rejected before any upload is spooled or any
filesystem operation takes place.

Legitimate names that merely *look* like internal names (for example
``.tmp.instrument-42``) are valid and must remain fully usable.
"""
from __future__ import annotations

import os
import re

# Deliberately permissive within one component: letters, digits, dot,
# underscore and hyphen.  No slash, no whitespace, no NUL.
_COMPONENT_RE = re.compile(r"[A-Za-z0-9_.-]{1,64}")
_RESERVED_COMPONENTS = frozenset({".", "..", ""})


def is_safe_component(name: object) -> bool:
    """True when ``name`` is a single ordinary path component."""
    if not isinstance(name, str) or not _COMPONENT_RE.fullmatch(name):
        return False
    if name in _RESERVED_COMPONENTS:
        return False
    # Defence in depth: the name must be its own basename and relative.
    return not os.path.isabs(name) and os.path.basename(name) == name


def is_valid_lane_name(name: object) -> bool:
    return is_safe_component(name)


def is_valid_sample_id(name: object) -> bool:
    return is_safe_component(name)
