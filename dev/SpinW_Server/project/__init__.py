# project/__init__.py


#################
#### imports ####
#################

import os

from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy
from flask_httpauth import HTTPBasicAuth

################
#### config ####
################

def _check_config_variables_are_set(config):

    assert config['SECRET_KEY'] is not None,\
           'SECRET_KEY is not set, set it in the production config file.'
    assert config['SECURITY_PASSWORD_SALT'] is not None,\
           'SECURITY_PASSWORD_SALT is not set, '\
           'set it in the production config file.'

    assert config['SQLALCHEMY_DATABASE_URI'] is not None,\
           'SQLALCHEMY_DATABASE_URI is not set, '\
           'set it in the production config file.'


app = Flask(__name__)


app.config.from_object(os.environ['APP_SETTINGS'])
config = app.config
_check_config_variables_are_set(config)

####################
#### extensions ####
####################

# extensions
auth = HTTPBasicAuth()
db = SQLAlchemy(app)

if config.get('USE_PYMATLAB'):
    import matlab.engine
    eng = matlab.engine.start_matlab()
else:
    from project.communication import tcip
    eng = tcip.tcip(config.get('DEPLOY_SERVER'), config.get('DEPLOY_PORT'))


####################
#### vairables  ####
####################
results = dict()
running = dict()

####################
#### blueprints ####
####################

from project.main.views import main_blueprint
from project.spinw.views import spinw_blueprint
from project.user.views import user_blueprint
app.register_blueprint(main_blueprint)
app.register_blueprint(user_blueprint)
app.register_blueprint(spinw_blueprint)
