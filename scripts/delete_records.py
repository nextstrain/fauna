import sys
import argparse
from rethinkdb import r
sys.path.append('')
from base.rethink_io import rethink_io
from vdb.download import rethinkdb_date_greater


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-db", "--database", help="database to delete from")
    parser.add_argument("-v", "--virus", default="flu", help="virus table to interact with")
    parser.add_argument("--filter", nargs="*", default=[], help="Filters for records to delete, i.e. inclusion_date: 2021-02-12")
    parser.add_argument("--match", nargs="*", default=[], help="Match for records to delete with regex pattern, e.g. accession:^EPIEPI")
    parser.add_argument('--interval', nargs="*", default=[], help="Select date interval of values for fields, e.g. assay_date:2019-09-03,2023-10-25")
    parser.add_argument("--preview", action="store_true", help="Preview records to be deleted without deleting from db.")

    args = parser.parse_args()

    rethink = rethink_io()
    rethink_host, auth_key = rethink.assign_rethink(rethink_host = None, auth_key = None)
    rethink.connect_rethink(args.database, rethink_host, auth_key)

    delete_filters = {}
    for delete_filter in args.filter:
        key, value = delete_filter.split(':')
        delete_filters[key] = value

    print(f"Delete filters: {delete_filters}")

    rethinkdb_command = r.table(args.virus).filter(delete_filters)

    delete_matches = {}
    for delete_match in args.match:
        field, pattern = delete_match.split(":")
        rethinkdb_command = rethinkdb_command.filter(lambda doc: doc[field].match(pattern))
        delete_matches[field] = pattern

    print(f"Delete matches: {delete_matches}")

    delete_intervals = {}
    for interval in args.interval:
        field, values = interval.split(':')
        older_date, newer_date = values.split(',')
        rethinkdb_command = rethinkdb_command.filter(lambda doc: rethinkdb_date_greater(doc[field].split('-'), older_date.split('-'), relaxed_interval=False))
        rethinkdb_command = rethinkdb_command.filter(lambda doc: rethinkdb_date_greater(newer_date.split('-'), doc[field].split('-'), relaxed_interval=False))
        delete_intervals[field] = (older_date, newer_date)

    print(f"Delete intervals: {delete_intervals}")

    if args.preview:
        filtered_records = rethinkdb_command.run()
        deletion_count = 0
        sources = set()
        for record in filtered_records:
            source = file_source = record.get('source', "")

            # The titer uploads using `elife_upload` appends a line counter to the end of the sources.
            if "_" in source:
                file_source, line = source.rsplit('_', 1)

            sources.add(file_source)
            deletion_count += 1

        print("Preview: selection would delete {} records".format(deletion_count))
        print("Sources of deleted records: {}".format(sources))
    else:
        deleted = rethinkdb_command.delete().run()
        print("Deleted {} records".format(deleted["deleted"]))
