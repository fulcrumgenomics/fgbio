'use strict';

// Tests for js/ncbi-assembly-report-to-sam-dict.js.  Run with: node --test _tests/
// (Jekyll does not publish directories that start with an underscore.)

const test = require("node:test");
const assert = require("node:assert/strict");
const C = require("../js/ncbi-assembly-report-to-sam-dict.js");

const Header = "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tGenBank-Accn\tRefSeq-Accn\tSequence-Length\tUCSC-style-name";

/** Builds a report from rows of [Sequence-Name, Sequence-Role, Assigned-Molecule, GenBank, RefSeq, length, UCSC]. */
function report(rows, { name = "Test1", eol = "\n" } = {}) {
  const lines = [`# Assembly name:  ${name}`, "# Organism name:  Testus testus", "#", "## not metadata: x", Header];
  return lines.concat(rows.map(r => r.join("\t"))).join(eol) + eol;
}

function build(text, options) {
  return C.buildSequenceDictionary(C.parseAssemblyReport(text), {
    primary: "UCSC-style-name",
    alternates: ["RefSeq-Accn", "GenBank-Accn", "Assigned-Molecule"],
    roles: [...C.DefaultSequenceRoles],
    sortByRole: true,
    assemblyId: "",
    includeReportTags: false,
    skipMissingAlternates: false,
    ...options
  });
}

const lines = result => result.text.trimEnd().split("\n");
const reasons = result => result.skipped.map(s => [s.name, s.reason]);

const Human = report([
  ["1", "assembled-molecule", "1", "CM000663.2", "NC_000001.11", "248956422", "chr1"],
  ["HSCHR1_CTG1_UNLOCALIZED", "unlocalized-scaffold", "1", "KI270706.1", "NT_187361.1", "175055", "chr1_KI270706v1_random"],
  ["HSCHR1_1_CTG3", "alt-scaffold", "1", "KI270762.1", "NT_187515.1", "354444", "chr1_KI270762v1_alt"],
  ["HG1_PATCH", "fix-patch", "1", "KN196472.1", "NW_009646194.1", "186494", "na"],
  ["HSCHRUN_RANDOM_CTG1", "unplaced-scaffold", "na", "KI270302.1", "NT_187395.1", "2274", "chrUn_KI270302v1"]
]);

test("parseAssemblyReport reads metadata, columns and records", () => {
  const r = C.parseAssemblyReport(Human);
  assert.equal(r.metadata.get("Assembly name"), "Test1");
  assert.equal(r.metadata.get("Organism name"), "Testus testus");
  assert.equal(r.metadata.has("not metadata"), false);
  assert.equal(r.columns.length, 7);
  assert.equal(r.records.length, 5);
  assert.equal(r.records[0]["UCSC-style-name"], "chr1");
});

test("parseAssemblyReport handles CRLF, a byte-order mark, and trailing tabs", () => {
  const text = "﻿" + report([["1", "assembled-molecule", "1", "A.1", "N.1", "10", "chr1\t\t"]], { eol: "\r\n" });
  const r = C.parseAssemblyReport(text);
  assert.equal(r.metadata.get("Assembly name"), "Test1");
  assert.equal(r.records[0]["UCSC-style-name"], "chr1");
});

test("parseAssemblyReport ignores comment lines after the records and a trailing tab on the header", () => {
  const text = `${Header}\t\n1\tassembled-molecule\t1\tA.1\tN.1\t10\tchr1\n# trailing comment\n#\tSequence-Role\n`;
  const r = C.parseAssemblyReport(text);
  assert.equal(r.columns.length, 7);
  assert.deepEqual(r.records, [{
    "Sequence-Name": "1", "Sequence-Role": "assembled-molecule", "Assigned-Molecule": "1", "GenBank-Accn": "A.1",
    "RefSeq-Accn": "N.1", "Sequence-Length": "10", "UCSC-style-name": "chr1"
  }]);
});

test("parseAssemblyReport rejects malformed reports", () => {
  assert.throws(() => C.parseAssemblyReport(""), /empty/);
  assert.throws(() => C.parseAssemblyReport("# foo: bar\n1\t2\n"), /Missing column header/);
  assert.throws(() => C.parseAssemblyReport("# foo: bar\n"), /Missing column header/);
  assert.throws(() => C.parseAssemblyReport(`${Header}\n`), /no sequences/);
  assert.throws(() => C.parseAssemblyReport("# Sequence-Name\tSequence-Role\n1\tassembled-molecule\n"), /missing the required column "Sequence-Length"/);
  const long = "x".repeat(500);
  const bad = `1\tassembled-molecule\t${long}`;
  assert.throws(() => C.parseAssemblyReport(`${Header}\n${bad}\n`),
    { message: `Line 2 has 3 column(s) but the header has 7: ${bad.slice(0, 100)}...` });
});

test("inferAssemblyId uses the Assembly name header, then the file name", () => {
  assert.equal(C.inferAssemblyId(C.parseAssemblyReport(Human).metadata, "x.txt"), "Test1");
  assert.equal(C.inferAssemblyId(new Map(), "GCF_000001405.40_GRCh38.p14_assembly_report.txt"), "GRCh38.p14");
  assert.equal(C.inferAssemblyId(new Map(), "/a/b/GRCh38.p13.assembly_report.txt"), "GRCh38.p13");
  assert.equal(C.inferAssemblyId(new Map(), ""), "");
});

test("rolesInReport lists default roles first, then other known roles, then unknown roles", () => {
  const r = C.parseAssemblyReport(report([
    ["a", "zeta-role", "na", "A.1", "na", "1", "na"],
    ["b", "fix-patch", "na", "B.1", "na", "1", "na"],
    ["c", "alt-scaffold", "na", "C.1", "na", "1", "na"],
    ["d", "assembled-molecule", "d", "D.1", "na", "1", "na"],
    ["e", "unplaced-scaffold", "na", "E.1", "na", "1", "na"]
  ]));
  assert.deepEqual(C.rolesInReport(r.records), ["assembled-molecule", "unplaced-scaffold", "alt-scaffold", "fix-patch", "zeta-role"]);
  assert.deepEqual([...C.roleCounts(r.records)], [["zeta-role", 1], ["fix-patch", 1], ["alt-scaffold", 1], ["assembled-molecule", 1], ["unplaced-scaffold", 1]]);
});

test("defaultColumns prefers UCSC names when they name every sequence in the default roles", () => {
  const r = C.parseAssemblyReport(Human);
  assert.deepEqual(C.defaultColumns(r.columns, r.records), {
    primary: "UCSC-style-name", alternates: ["RefSeq-Accn", "GenBank-Accn", "Assigned-Molecule"]
  });
});

test("defaultColumns falls back when UCSC names are incomplete (e.g. C. elegans, GRCm38 alt scaffolds)", () => {
  const worm = C.parseAssemblyReport(report([
    ["I", "assembled-molecule", "I", "BX284601.5", "NC_003279.8", "15072434", "na"],
    ["II", "assembled-molecule", "II", "BX284602.5", "NC_003280.10", "15279421", "na"],
    ["MtDNA", "assembled-molecule", "MtDNA", "X54252.1", "NC_001328.1", "13794", "chrM"]
  ]));
  assert.deepEqual(C.defaultColumns(worm.columns, worm.records), {
    primary: "RefSeq-Accn", alternates: ["GenBank-Accn", "Assigned-Molecule", "UCSC-style-name"]
  });
  const mouse = C.parseAssemblyReport(report([
    ["1", "assembled-molecule", "1", "CM000994.2", "NC_000067.6", "195471971", "chr1"],
    ["CHR_1_ALT", "alt-scaffold", "1", "KB469738.1", "NT_187001.1", "100", "na"],
    ["PATCH", "fix-patch", "1", "KB469739.1", "na", "100", "na"]
  ]));
  assert.equal(C.defaultColumns(mouse.columns, mouse.records).primary, "RefSeq-Accn");
  // No column names everything: use the one that names the most
  const partial = C.parseAssemblyReport(report([
    ["a", "unplaced-scaffold", "na", "A.1", "na", "1", "chrA"],
    ["b", "unplaced-scaffold", "na", "na", "na", "1", "chrB"],
    ["na", "unplaced-scaffold", "na", "C.1", "na", "1", "na"]
  ]));
  assert.equal(C.defaultColumns(partial.columns, partial.records).primary, "UCSC-style-name");
});

test("buildSequenceDictionary writes @HD, SN, LN, AS and AN, sorted by role", () => {
  const r = build(Human, { assemblyId: "Test1" });
  assert.deepEqual(lines(r), [
    "@HD\tVN:1.6",
    "@SQ\tSN:chr1\tLN:248956422\tAS:Test1\tAN:NC_000001.11,CM000663.2,1",
    "@SQ\tSN:chr1_KI270706v1_random\tLN:175055\tAS:Test1\tAN:NT_187361.1,KI270706.1",
    "@SQ\tSN:chrUn_KI270302v1\tLN:2274\tAS:Test1\tAN:NT_187395.1,KI270302.1",
    "@SQ\tSN:chr1_KI270762v1_alt\tLN:354444\tAS:Test1\tAN:NT_187515.1,KI270762.1"
  ]);
  assert.deepEqual(reasons(r), [["HG1_PATCH", 'sequence role "fix-patch" not selected']]);
  assert.equal(r.totalLength, 248956422 + 175055 + 2274 + 354444);
  assert.deepEqual(r.alternateColumns, ["RefSeq-Accn", "GenBank-Accn", "Assigned-Molecule"]);
});

test("buildSequenceDictionary follows the role order, or report order when not sorting", () => {
  const roles = ["alt-scaffold", "assembled-molecule"];
  assert.deepEqual(build(Human, { roles }).sequences.map(s => s.name), ["chr1_KI270762v1_alt", "chr1"]);
  assert.deepEqual(build(Human, { roles, sortByRole: false }).sequences.map(s => s.name), ["chr1", "chr1_KI270762v1_alt"]);
});

test("buildSequenceDictionary writes the report tags, omitting invalid values", () => {
  const r = build(Human, { roles: ["assembled-molecule"], includeReportTags: true });
  assert.deepEqual(lines(r)[1].split("\t"), [
    "@SQ", "SN:chr1", "LN:248956422", "AN:NC_000001.11,CM000663.2,1",
    "sr:assembled-molecule", "sn:1", "am:1", "ga:CM000663.2", "ra:NC_000001.11", "un:chr1"
  ]);
  const odd = build(report([["café", "assembled-molecule", "1", "A.1", "N.1", "10", "chr1"]]), { includeReportTags: true });
  assert.equal(lines(odd)[1], "@SQ\tSN:chr1\tLN:10\tAN:N.1,A.1,1\tsr:assembled-molecule\tam:1\tga:A.1\tra:N.1\tun:chr1");
  assert.deepEqual(odd.warnings, ['chr1: omitted tag sn because its value "caf\u00e9" is empty or not printable ASCII.']);
});

test("buildSequenceDictionary omits an AS that is not printable ASCII", () => {
  const r = build(Human, { roles: ["assembled-molecule"], assemblyId: "GRCh38–x" });
  assert.equal(lines(r)[1], "@SQ\tSN:chr1\tLN:248956422\tAN:NC_000001.11,CM000663.2,1");
  assert.deepEqual(r.warnings, ['The assembly identifier (AS) "GRCh38\u2013x" contains characters other than printable ASCII and was omitted.']);
});

test("buildSequenceDictionary only uses Assigned-Molecule for assembled molecules, as primary or alternate", () => {
  const r = build(Human, { primary: "Assigned-Molecule", alternates: ["UCSC-style-name"], roles: ["alt-scaffold", "assembled-molecule"] });
  assert.deepEqual(lines(r), ["@HD\tVN:1.6", "@SQ\tSN:1\tLN:248956422\tAN:chr1"]);
  assert.deepEqual(reasons(r).filter(([, reason]) => reason.includes("parent molecule")), [
    ["HSCHR1_1_CTG3", 'Assigned-Molecule names the parent molecule, not the sequence, for "alt-scaffold" sequences']
  ]);
});

test("buildSequenceDictionary reports an empty Assigned-Molecule primary as missing", () => {
  const r = build(report([["u", "unplaced-scaffold", "", "A.1", "na", "10", "chrU"]]), { primary: "Assigned-Molecule", alternates: [] });
  assert.deepEqual(reasons(r), [["u", 'missing a value for the primary column "Assigned-Molecule"']]);
});

test("buildSequenceDictionary skips invalid, duplicate, unnamed and bad-length sequences", () => {
  const r = build(report([
    ["a", "unplaced-scaffold", "na", "A.1", "na", "10", "chrA"],
    ["b", "unplaced-scaffold", "na", "B.1", "na", "10", "chrA"],
    ["c", "unplaced-scaffold", "na", "C.1", "na", "na", "chrC"],
    ["d", "unplaced-scaffold", "na", "D.1", "na", "0", "chrD"],
    ["e", "unplaced-scaffold", "na", "E.1", "na", "10", "na"],
    ["f", "unplaced-scaffold", "na", "F.1", "na", "10", "bad name"],
    ["g", "fix-patch", "na", "G.1", "na", "10", "na"]
  ]));
  assert.deepEqual(lines(r).slice(1), ["@SQ\tSN:chrA\tLN:10\tAN:A.1"]);
  assert.deepEqual(reasons(r), [
    ["g", 'sequence role "fix-patch" not selected'],
    ["chrA", "duplicate primary name"],
    ["chrC", 'invalid sequence length "na"'],
    ["chrD", 'invalid sequence length "0"'],
    ["e", 'missing a value for the primary column "UCSC-style-name"'],
    ["bad name", "not a valid SAM reference sequence name"]
  ]);
});

test("buildSequenceDictionary drops ambiguous alternate names from every sequence", () => {
  const r = build(report([
    ["a", "unplaced-scaffold", "na", "SHARED.1", "chrB", "10", "chrA"],
    ["b", "unplaced-scaffold", "na", "SHARED.1", "OK.1", "10", "chrB"],
    ["c", "unplaced-scaffold", "na", "C.1", "C.1", "10", "chrC"]
  ]));
  assert.deepEqual(lines(r).slice(1), [
    "@SQ\tSN:chrA\tLN:10",
    "@SQ\tSN:chrB\tLN:10\tAN:OK.1",
    "@SQ\tSN:chrC\tLN:10\tAN:C.1"
  ]);
  assert.deepEqual(r.warnings, [
    'chrA: dropped alternate name "chrB" from RefSeq-Accn because it is the primary name of another sequence.',
    'Dropped alternate name "SHARED.1" from chrA, chrB because it is shared by more than one sequence.'
  ]);
});

test("buildSequenceDictionary accepts alternate names allowed by the SAM spec's name rule", () => {
  const r = build(report([["a", "unplaced-scaffold", "na", "HLA-A*01:01:01:01", "na", "10", "chrA"]]));
  assert.equal(lines(r)[1], "@SQ\tSN:chrA\tLN:10\tAN:HLA-A*01:01:01:01");
  assert.equal(C.AlternateNameRegex.test("*bad"), false);
  assert.equal(C.AlternateNameRegex.test("a,b"), false);
});

test("buildSequenceDictionary frees the names of sequences skipped for having no alternates", () => {
  const r = build(report([
    ["c", "unplaced-scaffold", "na", "chrB", "C.1", "10", "chrC"],
    ["b", "unplaced-scaffold", "na", "na", "na", "10", "chrB"]
  ]), { skipMissingAlternates: true });
  assert.deepEqual(lines(r).slice(1), ["@SQ\tSN:chrC\tLN:10\tAN:C.1,chrB"]);
  assert.deepEqual(reasons(r), [["chrB", "no alternate names"]]);
  assert.deepEqual(r.warnings, []);
});

test("buildSequenceDictionary keeps a sequence whose only alternate is freed by a skipped sequence", () => {
  const r = build(report([
    ["c", "unplaced-scaffold", "na", "chrB", "na", "10", "chrC"],
    ["b", "unplaced-scaffold", "na", "na", "na", "10", "chrB"]
  ]), { skipMissingAlternates: true });
  assert.deepEqual(lines(r).slice(1), ["@SQ\tSN:chrC\tLN:10\tAN:chrB"]);
  assert.deepEqual(reasons(r), [["chrB", "no alternate names"]]);
  assert.deepEqual(r.warnings, []);
});

test("buildSequenceDictionary skips one of two sequences whose only alternates are each other's names", () => {
  const r = build(report([
    ["c", "unplaced-scaffold", "na", "chrB", "na", "10", "chrC"],
    ["b", "unplaced-scaffold", "na", "chrC", "na", "10", "chrB"]
  ]), { skipMissingAlternates: true });
  assert.deepEqual(lines(r).slice(1), ["@SQ\tSN:chrB\tLN:10\tAN:chrC"]);
  assert.deepEqual(reasons(r), [["chrC", "no alternate names"]]);
});

test("buildSequenceDictionary keeps a sequence that waits on a cycle once the cycle is skipped", () => {
  const r = build(report([
    ["a", "unplaced-scaffold", "na", "chrB", "na", "10", "chrA"],
    ["b", "unplaced-scaffold", "na", "chrA", "na", "10", "chrB"],
    ["c", "unplaced-scaffold", "na", "chrA", "na", "10", "chrC"]
  ]), { skipMissingAlternates: true });
  assert.deepEqual(lines(r).slice(1), ["@SQ\tSN:chrC\tLN:10\tAN:chrA"]);
  assert.deepEqual(reasons(r), [["chrA", "no alternate names"], ["chrB", "no alternate names"]]);
});

test("buildSequenceDictionary breaks cycles one sequence at a time", () => {
  const r = build(report([
    ["a", "unplaced-scaffold", "na", "chrB", "na", "10", "chrA"],
    ["b", "unplaced-scaffold", "na", "chrA", "chrD", "10", "chrB"],
    ["d", "unplaced-scaffold", "na", "chrE", "na", "10", "chrD"],
    ["e", "unplaced-scaffold", "na", "chrD", "na", "10", "chrE"]
  ]), { skipMissingAlternates: true });
  assert.deepEqual(lines(r).slice(1), ["@SQ\tSN:chrB\tLN:10\tAN:chrA", "@SQ\tSN:chrD\tLN:10\tAN:chrE"]);
  assert.deepEqual(reasons(r), [["chrA", "no alternate names"], ["chrE", "no alternate names"]]);
});

test("buildSequenceDictionary does not wait on a name another sequence with alternates also claims", () => {
  const r = build(report([
    ["a", "unplaced-scaffold", "na", "chrB", "na", "10", "chrA"],
    ["b", "unplaced-scaffold", "na", "chrA", "na", "10", "chrB"],
    ["c", "unplaced-scaffold", "na", "chrB", "C.1", "10", "chrC"]
  ]), { skipMissingAlternates: true });
  assert.deepEqual(lines(r).slice(1), ["@SQ\tSN:chrB\tLN:10\tAN:chrA", "@SQ\tSN:chrC\tLN:10\tAN:C.1"]);
  assert.deepEqual(reasons(r), [["chrA", "no alternate names"]]);
});

test("buildSequenceDictionary gives a shared alternate to the remaining claimant after skipping the other", () => {
  const r = build(report([
    ["c", "unplaced-scaffold", "na", "X.1", "na", "10", "chrC"],
    ["d", "unplaced-scaffold", "na", "X.1", "D.1", "10", "chrD"]
  ]), { skipMissingAlternates: true });
  assert.deepEqual(lines(r).slice(1), ["@SQ\tSN:chrD\tLN:10\tAN:D.1,X.1"]);
  assert.deepEqual(reasons(r), [["chrC", "no alternate names"]]);
  assert.deepEqual(r.warnings, []);
});

test("parseAccession, accessionDirectoryUrl, findAssemblyDirectory and assemblyReportUrl", () => {
  assert.deepEqual(C.parseAccession(" gcf_000001405.40 "), { prefix: "GCF", digits: "000001405", version: 40 });
  assert.deepEqual(C.parseAccession("GCA_000001635"), { prefix: "GCA", digits: "000001635", version: null });
  assert.equal(C.parseAccession("GRCh38"), null);
  const acc = C.parseAccession("GCF_000001635");
  assert.equal(C.accessionDirectoryUrl(acc), "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/");
  const hrefs = ["/genomes/all/GCF/000/001/", "GCF_000001635.26_GRCm38.p6/", "GCF_000001635.27_GRCm39/", "GCF_000001635.9_X/", "GCA_000001635.30_Y/"];
  assert.equal(C.findAssemblyDirectory(hrefs, acc), "GCF_000001635.27_GRCm39");
  assert.equal(C.findAssemblyDirectory(hrefs, C.parseAccession("GCF_000001635.26")), "GCF_000001635.26_GRCm38.p6");
  assert.equal(C.findAssemblyDirectory(hrefs, C.parseAccession("GCF_000001635.2")), null);
  assert.equal(C.assemblyReportUrl(acc, "GCF_000001635.27_GRCm39"),
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt");
});

test("withThousands", () => {
  assert.equal(C.withThousands(3209286105), "3,209,286,105");
  assert.equal(C.withThousands(12), "12");
});
