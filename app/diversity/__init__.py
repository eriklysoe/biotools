"""Diversity blueprint – Alpha & Beta diversity analysis tool."""

import os

from flask import Blueprint, request, Response

diversity_bp = Blueprint("diversity", __name__)

ADMIN_USER = os.environ.get("ADMIN_USER", "")
ADMIN_PASS = os.environ.get("ADMIN_PASS", "")


@diversity_bp.before_request
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
