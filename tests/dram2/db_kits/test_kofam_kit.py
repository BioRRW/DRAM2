


def test_kofam(runner):
  result = runner.invoke(hello, ['Peter'])
  assert result.exit_code == 0
  assert result.output == 'Hello Peter!\n'

def test_kofam_hmmscan_formater_dbcan():
    output_expt = pd.DataFrame(
        {
            "bin_1.scaffold_1": ["K00001", "KO1; description 1"],
            "bin_1.scaffold_2": ["K00002", "KO2; description 2"],
        },
        index=["ko_id", "kegg_hit"],
    ).T
    input_b6 = os.path.join("dram2", "annotate", "tests", "data", "unformatted_kofam.b6")
    hits = parse_hmmsearch_domtblout(input_b6)
    output_rcvd = kofam_hmmscan_formater(
        hits=hits,
        hmm_info_path=os.path.join("dram2", "annotate", "tests", "data", "hmm_thresholds.txt"),
        top_hit=True,
        use_dbcan2_thresholds=True,
    )
    output_rcvd.sort_index(inplace=True)
    output_expt.sort_index(inplace=True)
    assert output_rcvd.equals(
        output_expt
    ), "Error in kofam_hmmscan_formater with dbcam cuts"


def test_kofam_hmmscan_formater():
    output_expt = pd.DataFrame(
        {"bin_1.scaffold_2": ["K00002", "KO2; description 2"]},
        index=["ko_id", "kegg_hit"],
    ).T
    input_b6 = os.path.join("dram2", "annotate", "tests", "data", "unformatted_kofam.b6")
    hits = parse_hmmsearch_domtblout(input_b6)
    output_rcvd = kofam_hmmscan_formater(
        hits=hits,
        hmm_info_path=os.path.join("dram2", "annotate", "tests", "data", "hmm_thresholds.txt"),
        top_hit=True,
        use_dbcan2_thresholds=False,
    )
    output_rcvd.sort_index(inplace=True)
    output_expt.sort_index(inplace=True)
    assert output_rcvd.equals(output_expt), "Error in kofam_hmmscan_formater"

