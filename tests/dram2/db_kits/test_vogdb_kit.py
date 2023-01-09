


def test_vogdb_hmmscan_formater():
    # TODO Not Done
    output_expt = pd.DataFrame(
        {"bin_1.scaffold_1": "VOG00001", "bin_1.scaffold_2": "VOG00002"},
        index=["vogdb_id"],
    ).T
    input_b6 = os.path.join("dram2", "annotate", "tests", "data", "unformatted_vogdb.b6")
    hits = parse_hmmsearch_domtblout(input_b6)
    output_rcvd = vogdb_hmmscan_formater(hits=hits, db_name="vogdb")
    output_rcvd.sort_index(inplace=True)
    output_expt.sort_index(inplace=True)
    assert output_rcvd.equals(output_expt), "Error in vogdb_hmmscan_formater"
