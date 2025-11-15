import sys
import os
import logging

logging.basicConfig(stream=sys.stderr)

# Make sure the app package is on the path
sys.path.insert(0, "/var/www/FlaskApp")
sys.path.insert(0, "/var/www/FlaskApp/FlaskApp")

# Import the Flask app object
from FlaskApp import app as application

# Optional: set secret key (better to use env var in real deployments)
application.secret_key = os.environ.get("FLASK_SECRET_KEY", "change-me")
