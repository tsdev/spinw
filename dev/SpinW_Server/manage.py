# manage.py


import os

# from flask_script import Manager
# from flask_migrate import Migrate, MigrateCommand
#
# COV = coverage.coverage(
#         branch=True,
#         include='project/*',
#         omit=['*/__init__.py', '*/config/*']
#     )
# COV.start()

from project import app, db
from project.models import User, UserJobs

# migrate = Migrate(app, db)
# manager = Manager(app)

# migrations
# manager.add_command('db', MigrateCommand)

#
# @manager.command
# def test():
#     """Runs the unit tests without coverage."""
#     tests = unittest.TestLoader().discover('tests')
#     result = unittest.TextTestRunner(verbosity=2).run(tests)
#     if result.wasSuccessful():
#         return 0
#     else:
#         return 1
#
#
# @manager.command
# def cov():
#     """Runs the unit tests with coverage."""
#     tests = unittest.TestLoader().discover('tests')
#     unittest.TextTestRunner(verbosity=2).run(tests)
#     COV.stop()
#     COV.save()
#     print('Coverage Summary:')
#     COV.report()
#     basedir = os.path.abspath(os.path.dirname(__file__))
#     covdir = os.path.join(basedir, 'tmp/coverage')
#     COV.html_report(directory=covdir)
#     print('HTML version: file://%s/index.html' % covdir)
#     COV.erase()
#
#
# @manager.command
# def create_db():
#     """Creates the db tables."""
#     db.create_all()
#
#
# @manager.command
# def drop_db():
#     """Drops the db tables."""
#     db.drop_all()
#
#
# @manager.command
# def create_admin():
#     """Creates the admin user."""
#     db.session.add(User(
#         email="ad@min.com",
#         password="admin",
#         admin=True,
#         confirmed=True,
#         confirmed_on=datetime.datetime.now())
#     )
#     db.session.commit()


if __name__ == '__main__':
    # manager.run()
    app.config.from_object(os.environ['APP_SETTINGS'])
    if app.config['SQLALCHEMY_DATABASE_URI'] is 'sqlite:///:memory:':
        db.create_all()
    else:
        if not os.path.exists(app.config['SQLALCHEMY_DATABASE_URI'].split('///')[1]):
            db.create_all()
    UserJobs.__table__.drop(db.engine)
    UserJobs.__table__.create(db.engine)

    ####################
    #### threading  ####
    ####################

    from status import Status

    status_thread = Status(2)
    status_thread.start()
    app.run(debug=True)