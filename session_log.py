import os
import sys
import datetime

LOG_PATH = os.path.join(os.path.dirname(__file__),
                        'docs', 'decisions', 'parametricviz_session_log.md')

sys.path.insert(0, os.path.expanduser('~/GitHub_Projects/shared'))
from project_status_updater import update_project_status


def log_session(title, description, next_step,
                current_file, open_tasks, next_session, notes=''):
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

    entry = (
        f'\n---\n\n'
        f'**{timestamp} — {title}**\n\n'
        f'{description}\n\n'
        f'**Next:** {next_step}\n'
    )

    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(entry)

    print(f'  session log updated — {title}')

    update_project_status(
        project_id   = 'parametric_viz',
        display_name = 'Parametric Viz',
        current_file = current_file,
        open_tasks   = open_tasks,
        next_session = next_session,
        notes        = notes,
    )
