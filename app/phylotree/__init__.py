"""PhyloTree blueprint – phylogenetic tree builder."""

import os
from functools import wraps

from flask import Blueprint, request, Response

phylotree_bp = Blueprint("phylotree", __name__)

ADMIN_USER = os.environ.get("ADMIN_USER", "")
ADMIN_PASS = os.environ.get("ADMIN_PASS", "")


@phylotree_bp.before_request
def _check_auth():
    if not ADMIN_USER or not ADMIN_PASS:
        return
    auth = request.authorization
    if not auth or auth.username != ADMIN_USER or auth.password != ADMIN_PASS:
        return Response(
            "Login required.", 401,
            {"WWW-Authenticate": 'Basic realm="BioTools"'},
        )


from . import routes  # noqa: E402, F401
