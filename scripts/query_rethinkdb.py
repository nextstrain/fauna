import sys
import argparse
from rethinkdb import r
sys.path.append('')
from base.rethink_io import rethink_io


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-db", "--database", help="database to delete from")
    parser.add_argument("-v", "--virus", default="flu", help="virus table to interact with")
    parser.add_argument("--filter", nargs="*", default=[], help="Filters for records, e.g inclusion_date: 2021-02-12")

    args = parser.parse_args()

    rethink = rethink_io()
    rethink_host, auth_key = rethink.assign_rethink(rethink_host = None, auth_key = None)
    rethink.connect_rethink(args.database, rethink_host, auth_key)

    filters = {}
    for delete_filter in args.filter:
        key, value = delete_filter.split(':')
        filters[key] = value

    print(f"Filters: {filters}")

    rethinkdb_command = r.table(args.virus).filter(filters)

    filtered_records = rethinkdb_command.run()

    for record in filtered_records:
        print(record)
