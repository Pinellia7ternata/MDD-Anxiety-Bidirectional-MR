"""Path resolution so the revision scripts run from the repository root."""
import os

_INPUTS = (os.environ.get('MR_DATA_DIR', ''), 'data', '../data',
           'JAD_Submission_Package_MDD_Anxiety_v2/06_Source_Data_CSV',
           '_TRASH_20260731/过程文档归档_20260804/06_上游源数据与工作目录/MR_results')


def locate(name):
    for base in _INPUTS:
        if base:
            p = os.path.join(base, name)
            if os.path.exists(p):
                return p
    raise FileNotFoundError(f'{name} not found in {_INPUTS}')


def output(name):
    out = os.environ.get('MR_OUT_DIR') or ('.' if os.path.isdir('data') else
                                           'JAD_Submission_Package_MDD_Anxiety_v2/08_Revision_20261102')
    os.makedirs(out, exist_ok=True)
    return os.path.join(out, name)
