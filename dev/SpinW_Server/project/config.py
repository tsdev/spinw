# project/config.py

import os
import configparser

basedir = os.path.abspath(os.path.dirname(__file__))


def _get_bool_env_var(varname, default=None):

    value = os.environ.get(varname, default)

    if value is None:
        return False
    elif isinstance(value, str) and value.lower() == 'false':
        return False
    elif bool(value) is False:
        return False
    else:
        return bool(value)


class BaseConfig(object):
    """Base configuration."""

    # main config
    SECRET_KEY = 'the quick brown fox jumps over the lazy dog'
    SECURITY_PASSWORD_SALT = 'blah blah blah'
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    UPLOAD_FOLDER = 'Uploaded'
    LDAP_SERVER = 'd.psi.ch'
    SERVER_VERSION = '2.0.0'
    USE_LDAP = False
    USE_PYMATLAB = True

class DevelopmentConfig(BaseConfig):
    """Development configuration."""
    DEPLOY_SERVER= '127.0.0.1'
    DEPLOY_PORT = 13001
    SQLALCHEMY_DATABASE_URI = 'sqlite:///' + os.path.join(basedir, 'db.sqlite')
    TOKEN_DURATION = 600  # This is 10 minutes


class TestingConfig(BaseConfig):
    """Testing configuration."""
    DEPLOY_SERVER= '127.0.0.1'
    DEPLOY_PORT = 13001
    SQLALCHEMY_DATABASE_URI = 'sqlite://:memory:'
    TOKEN_DURATION = 300  # This is 5 minutes


class ProductionConfig(BaseConfig):
    """Production configuration."""
    SECRET_KEY = None
    SECURITY_PASSWORD_SALT = None

    SQLALCHEMY_DATABASE_URI = None

    # production config takes precedence over env variables

    # production config file at ./project/config/production.cfg
    config_path = os.path.join(basedir, 'config', 'production.cfg')

    # if config file exists, read it:
    if os.path.isfile(config_path):
        config = configparser.ConfigParser()

        with open(config_path) as configfile:
            config.readfp(configfile)

        SECRET_KEY = config.get('keys', 'SECRET_KEY')
        SECURITY_PASSWORD_SALT = config.get('keys', 'SECRET_KEY')

        UPLOAD_FOLDER = config.get('server', 'UPLOAD_FOLDER')
        DEPLOY_SERVER = config.get('server', 'DEPLOY_SERVER')
        DEPLOY_PORT = config.get('server', 'DEPLOY_PORT')
        TOKEN_DURATION = config.get('server', 'TOKEN_DURATION')
        LDAP_SERVER = config.get('server', 'LDAP_SERVER')

        # database URI
        SQLALCHEMY_DATABASE_URI = config.get('db', 'SQLALCHEMY_DATABASE_URI')

        # Connections
        USE_LDAP = config.get('conn','USE_LDAP')
        USE_PYMATLAB = config.get('conn','USE_PYMATLAB')

