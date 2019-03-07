from models import engine, paper_status, paper_status_cnt

with engine.connect() as conn:
    conn.execute('truncate var_pmid')
