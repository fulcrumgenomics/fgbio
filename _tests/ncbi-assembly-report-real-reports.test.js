'use strict';

// Tests js/ncbi-assembly-report-to-sam-dict.js against real NCBI assembly reports in _tests/data.  Run with:
// node --test _tests/
//
// GRCh38.p13.CollectAlternateContigNames.dict was made with fgbio 4.1.0:
//   fgbio CollectAlternateContigNames -i GRCh38.p13.assembly_report.txt -o GRCh38.p13.CollectAlternateContigNames.dict \
//     -p UcscName -a RefSeqAccession GenBankAccession AssignedMolecule \
//     -s AssembledMolecule UnlocalizedScaffold UnplacedScaffold AltScaffold \
//     --sort-by-sequencing-role true --skip-missing-alternates false

const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const C = require("../js/ncbi-assembly-report-to-sam-dict.js");

const DataDir = path.join(__dirname, "data");
const read = name => fs.readFileSync(path.join(DataDir, name), "utf8");

/** Sorts the tags after the record type so lines compare regardless of tag order. */
const canonical = line => {
  const [type, ...fields] = line.split("\t");
  return [type, ...fields.sort()].join("\t");
};

test("GRCh38.p13 matches fgbio CollectAlternateContigNames with the same settings", () => {
  const report = C.parseAssemblyReport(read("GRCh38.p13.assembly_report.txt"));
  const result = C.buildSequenceDictionary(report, {
    primary: "UCSC-style-name",
    alternates: ["RefSeq-Accn", "GenBank-Accn", "Assigned-Molecule"],
    roles: [...C.DefaultSequenceRoles],
    sortByRole: true,
    includeReportTags: true
  });
  const expected = read("GRCh38.p13.CollectAlternateContigNames.dict").trim().split("\n").map(canonical);
  assert.deepEqual(result.text.trim().split("\n").map(canonical), expected);
  assert.equal(result.sequences.length, 455);
});

for (const file of fs.readdirSync(DataDir).filter(f => f.startsWith("GCF_") && f.endsWith("_assembly_report.txt")).sort()) {
  test(`${file}: the page's default settings keep every sequence in the default roles, with unique names`, () => {
    const report = C.parseAssemblyReport(read(file));
    const columns = C.defaultColumns(report.columns, report.records);
    const present = C.rolesInReport(report.records);
    let roles = present.filter(role => C.DefaultSequenceRoles.includes(role));
    if (roles.length === 0) roles = present;
    const result = C.buildSequenceDictionary(report, {
      ...columns,
      roles,
      sortByRole: true,
      assemblyId: C.inferAssemblyId(report.metadata, file),
      includeReportTags: true,
      skipMissingAlternates: false
    });
    const inRoles = report.records.filter(record => roles.includes(record["Sequence-Role"])).length;
    assert.equal(result.sequences.length, inRoles);
    assert.deepEqual(result.warnings, []);

    const names = new Set();
    for (const line of result.text.trim().split("\n").slice(1)) {
      const tags = new Map(line.split("\t").slice(1).map(field => [field.slice(0, 2), field.slice(3)]));
      for (const name of [tags.get("SN"), ...(tags.get("AN") || "").split(",").filter(Boolean)]) {
        assert.ok(!names.has(name), `duplicate name ${name}`);
        names.add(name);
      }
    }
  });
}
