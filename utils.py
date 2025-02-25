import json
from calendar import monthrange
from datetime import datetime
import hashlib

def increment_months(start_date, months):
    # Compute the total number of months since year 0
    total_months = start_date.year * 12 + (start_date.month - 1) + months
    new_year = total_months // 12
    new_month = total_months % 12 + 1

    # Handle end-of-month cases by adjusting the day if necessary
    last_day_of_new_month = monthrange(new_year, new_month)[1]
    new_day = min(start_date.day, last_day_of_new_month)

    return datetime(new_year, new_month, new_day)

def generate_hash(obj):
    obj_string = json.dumps(obj).encode('utf8')
    obj_hash = hashlib.sha1(obj_string).hexdigest()[:12] # truncated for shorter hash
    print('generating hash from object:', obj_string, '->', obj_hash)
    return obj_hash
