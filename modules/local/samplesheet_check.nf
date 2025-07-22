process SAMPLESHEET_CHECK {
    tag "${samplesheet}"
    label 'process_single'
    container "ghcr.io/btc-absolute:latest"

    input:
    path samplesheet
    val phase

    output:
    path 'samplesheet_utf8.csv'    , emit: samplesheet_utf8

    script: 
    """
    python - <<EOF
    import csv
    import sys

    REQUIRED_COLUMNS = {
        '1': ['sample', 'seg_path', 'indel_path', 'snp_path'],
        '2': ['sample', 'seg_path', 'indel_path', 'snp_path', 'purity', 'ploidy', 'rdata_path']
    }

    def samplesheet(input_path, phase):
        required_fields = REQUIRED_COLUMNS.get(phase)
        if not required_fields:
            print(f"Invalid phase: {phase}. Must be 1 or 2.")
            sys.exit(1)
        try:
            with open(input_path, newline='', encoding='utf-8') as infile:
                reader = csv.DictReader(infile)
                fieldnames = reader.fieldnames

                # Validate required columns
                missing_cols = [col for col in required_fields if col not in fieldnames]
                if missing_cols:
                    print(f"ERROR: Samplesheet is missing required columns for Phase {phase}: {missing_cols}")
                    sys.exit(1)

                rows = list(reader)

            # Warn about rows missing values in required fields
            for i, row in enumerate(rows):
                missing_vals = [col for col in required_fields if not row[col].strip()]
                if missing_vals:
                    print(f"WARNING: Row {i + 2} is missing values in required fields: {missing_vals}")

            # Write UTF-8 BOM version for Excel compatibility
            with open('samplesheet_utf8.csv', 'w', newline='', encoding='utf-8-sig') as outfile:
                writer = csv.DictWriter(outfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)

        except Exception as e:
            print(f"ERROR: Exception during parsing: {e}", file=sys.stderr)
            sys.exit(1)

    if __name__ == "__main__":
        samplesheet("${samplesheet}", "${phase}")
    EOF
    """

    stub:
    """
    #!/bin/bash

    touch samplesheet_utf8.csv
    """
}
