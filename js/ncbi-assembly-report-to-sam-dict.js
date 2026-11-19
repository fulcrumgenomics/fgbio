'use strict';

/**
 * Builds a SAM sequence dictionary from an NCBI assembly report.  Used by ncbi-assembly-report-to-sam-dict.html.
 *
 * This file has no DOM access.  In Node it exports its functions (see _tests/ncbi-assembly-report-to-sam-dict.test.js);
 * in a browser the top-level declarations are shared with the page's own script.
 */

/** The value used in assembly reports for a missing column value. */
const MissingColumnValue = "na";

/**
 * Columns in an NCBI assembly report that can be used as sequence names, along with the lowercase SAM tag used to
 * store the column's value (the same tags as fgbio's CollectAlternateContigNames).
 */
const NameColumns = Object.freeze([
  Object.freeze({ key: "Sequence-Name",     tag: "sn" }),
  Object.freeze({ key: "Assigned-Molecule", tag: "am" }),
  Object.freeze({ key: "GenBank-Accn",      tag: "ga" }),
  Object.freeze({ key: "RefSeq-Accn",       tag: "ra" }),
  Object.freeze({ key: "UCSC-style-name",   tag: "un" })
]);
const SequenceRoleColumn   = Object.freeze({ key: "Sequence-Role", tag: "sr" });
const SequenceLengthKey    = "Sequence-Length";
const AssignedMoleculeKey  = "Assigned-Molecule";
const AssembledMoleculeRole = "assembled-molecule";

/** The order in which name columns are preferred for the default primary name (`SN`). */
const PrimaryColumnPreference = Object.freeze(["UCSC-style-name", "RefSeq-Accn", "GenBank-Accn", "Sequence-Name", "Assigned-Molecule"]);

/** The order in which alternate name columns are offered, and written to `AN`. */
const AlternateColumnOrder = Object.freeze(["RefSeq-Accn", "GenBank-Accn", "Assigned-Molecule", "UCSC-style-name", "Sequence-Name"]);

/** Alternate name columns selected by default (when present and not the primary). */
const DefaultAlternateColumns = Object.freeze(["RefSeq-Accn", "GenBank-Accn", "Assigned-Molecule", "UCSC-style-name"]);

/** Roles included by default. */
const DefaultSequenceRoles = Object.freeze([
  "assembled-molecule", "unlocalized-scaffold", "unplaced-scaffold", "alt-scaffold"
]);

/** Known sequence roles in the order they are listed: the default roles first, then the remaining known roles. */
const KnownSequenceRoles = Object.freeze(DefaultSequenceRoles.concat(["fix-patch", "novel-patch"]));

/**
 * Valid reference sequence name (`@SQ SN`), from the SAM specification.  The specification and htsjdk apply the same
 * rule to each alternate name in `@SQ AN`; commas are not allowed, so joining alternate names with ',' is unambiguous.
 */
const SequenceNameRegex = /^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$/;
const AlternateNameRegex = SequenceNameRegex;
/** Valid header tag value (e.g. `AS`, or a lowercase report tag): printable ASCII, including spaces. */
const TagValueRegex = /^[ -~]+$/;
/** The maximum reference sequence length allowed by the SAM specification. */
const MaxSequenceLength = 2147483647;
/** The maximum number of characters of an offending line to echo in an error message. */
const MaxEchoedLineLength = 100;

const NumberFormat = new Intl.NumberFormat("en-US");

/** Returns the name column with the given key, or undefined. */
function nameColumn(key) {
  return NameColumns.find(c => c.key === key);
}

/** Truncates text to at most `max` characters, adding an ellipsis when truncated. */
function truncate(text, max) {
  return text.length <= max ? text : `${text.slice(0, max)}...`;
}

/**
 * Parses the text of an NCBI assembly report.
 *
 * Header lines start with '#'.  Lines of the form `# Key: value` are collected into `metadata`.  The column header is
 * the '#' line that names the columns (it starts with `# Sequence-Name`).  Extra empty trailing fields (from trailing
 * tabs) are ignored.
 *
 * Returns {metadata: Map<string,string>, columns: string[], records: Array<Object>}, where each record maps a column
 * name to its value.  Throws an Error with a descriptive message if the report is malformed.
 */
function parseAssemblyReport(text) {
  if (typeof text !== "string" || text.trim() === "") {
    throw new Error("The assembly report is empty.");
  }
  const lines = text.replace(/^﻿/, "").split(/\r?\n|\r/);
  const metadata = new Map();
  let columns = null;
  const records = [];

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    if (line.trim() === "") continue;
    if (line.startsWith("#")) {
      if (records.length > 0) continue; // trailing comment lines
      const body = line.replace(/^#+\s*/, "");
      if (body.startsWith(`${NameColumns[0].key}\t`) || body.split("\t").includes(SequenceRoleColumn.key)) {
        columns = body.split("\t").map(c => c.trim());
        while (columns.length > 0 && columns[columns.length - 1] === "") columns.pop();
      } else if (!line.startsWith("##")) {
        const m = body.match(/^([^:\t]+):\s*(.*?)\s*$/);
        if (m && !metadata.has(m[1].trim())) metadata.set(m[1].trim(), m[2]);
      }
      continue;
    }
    if (columns === null) {
      throw new Error(`Missing column header line (expected a line starting with "# ${NameColumns[0].key}") before line ${i + 1}.`);
    }
    const fields = line.split("\t").map(f => f.trim());
    while (fields.length > columns.length && fields[fields.length - 1] === "") fields.pop();
    if (fields.length !== columns.length) {
      throw new Error(`Line ${i + 1} has ${fields.length} column(s) but the header has ${columns.length}: ${truncate(line, MaxEchoedLineLength)}`);
    }
    const record = {};
    for (let j = 0; j < columns.length; j++) record[columns[j]] = fields[j];
    records.push(record);
  }

  if (columns === null) {
    throw new Error(`Missing column header line (expected a line starting with "# ${NameColumns[0].key}").`);
  }
  for (const required of [SequenceRoleColumn.key, SequenceLengthKey]) {
    if (!columns.includes(required)) throw new Error(`The column header is missing the required column "${required}".`);
  }
  if (!NameColumns.some(c => columns.includes(c.key))) {
    throw new Error(`The column header has none of the name columns: ${NameColumns.map(c => c.key).join(", ")}.`);
  }
  if (records.length === 0) {
    throw new Error("The assembly report has no sequences.");
  }
  return { metadata, columns, records };
}

/**
 * Infers the genome assembly identifier (`AS`) from the report's `Assembly name` header, falling back to the file
 * name (e.g. `GCF_000001405.40_GRCh38.p14_assembly_report.txt` yields `GRCh38.p14`).  Returns "" if neither is known.
 */
function inferAssemblyId(metadata, fileName) {
  const name = metadata && metadata.get("Assembly name");
  if (name && name.trim() !== "") return name.trim();
  if (fileName) {
    let base = String(fileName).split(/[\\/]/).pop();
    base = base.replace(/\.(txt|tsv)$/i, "").replace(/[._]assembly[._]report$/i, "");
    base = base.replace(/^GC[AF]_\d{9}\.\d+_/, "");
    return base;
  }
  return "";
}

/** Returns the distinct sequence roles in the report: known roles in `KnownSequenceRoles` order, then unknown roles sorted. */
function rolesInReport(records) {
  const present = new Set(records.map(r => r[SequenceRoleColumn.key]));
  const known = KnownSequenceRoles.filter(r => present.has(r));
  const unknown = [...present].filter(r => !KnownSequenceRoles.includes(r)).sort();
  return known.concat(unknown);
}

/** Returns a Map from sequence role to the number of records with that role. */
function roleCounts(records) {
  const counts = new Map();
  for (const r of records) {
    const role = r[SequenceRoleColumn.key];
    counts.set(role, (counts.get(role) || 0) + 1);
  }
  return counts;
}

/**
 * Returns the record's name in the given column, or null if it has none.  Assigned-Molecule only names
 * assembled-molecule sequences: for other roles it names the parent molecule (e.g. chr1_KI270706v1_random -> 1), so
 * it is treated as missing.
 */
function nameValue(record, column) {
  const value = record[column];
  if (value === undefined || value === "" || value === MissingColumnValue) return null;
  if (column === AssignedMoleculeKey && record[SequenceRoleColumn.key] !== AssembledMoleculeRole) return null;
  return value;
}

/** Returns a label for a record in messages: its primary name if it has one, otherwise its Sequence-Name. */
function recordLabel(record, primary) {
  return nameValue(record, primary) || nameValue(record, NameColumns[0].key) || "(unnamed)";
}

/** Returns the name columns present in the report that have at least one name. */
function usableNameColumns(columns, records) {
  return NameColumns.filter(c => columns.includes(c.key) && records.some(r => nameValue(r, c.key) !== null));
}

/**
 * Returns the default primary and alternate columns.
 *
 * The primary is the first column in `PrimaryColumnPreference` that names every sequence in the roles selected by
 * default (all roles if none of the default roles are present); if no column names them all, the column that names
 * the most is used (ties go to the preferred column).  The alternates are the `DefaultAlternateColumns` that have at
 * least one name, other than the primary.
 */
function defaultColumns(columns, records) {
  const presentDefaults = new Set(DefaultSequenceRoles.filter(r => records.some(rec => rec[SequenceRoleColumn.key] === r)));
  const selected = presentDefaults.size > 0 ? records.filter(r => presentDefaults.has(r[SequenceRoleColumn.key])) : records;
  let primary = null;
  let bestCoverage = -1;
  for (const key of PrimaryColumnPreference.filter(k => columns.includes(k))) {
    const coverage = selected.filter(r => nameValue(r, key) !== null).length;
    if (coverage > bestCoverage) { primary = key; bestCoverage = coverage; }
    if (coverage === selected.length) break;
  }
  const usable = usableNameColumns(columns, records).map(c => c.key);
  const alternates = DefaultAlternateColumns.filter(k => k !== primary && usable.includes(k));
  return { primary, alternates };
}

/**
 * Builds the sequence dictionary.
 *
 * options:
 *   primary:               the column for the primary name (`SN`)
 *   alternates:            the columns for the alternate names (`AN`), in order (the primary is ignored if present)
 *   roles:                 the sequence roles to include, in sort order
 *   sortByRole:            sort sequences by the order of `roles` (stable), otherwise keep report order
 *   assemblyId:            the `AS` value; omitted when empty
 *   includeReportTags:     add `sr`, `sn`, `am`, `ga`, `ra`, `un` tags with the assembly report values
 *   skipMissingAlternates: skip sequences with no alternate names
 *
 * Every name in the output is unique across all `SN` and `AN` values.  An alternate name that is another output
 * sequence's primary name, or that is shared by more than one output sequence, is dropped (from all of them) with a
 * warning, since it would be ambiguous.
 *
 * Returns {sequences, skipped, warnings, text, totalLength, alternateColumns}, where each sequence is
 * {name, length, role, alternates: string[], altByColumn: Map<column, name>, record} and each skipped entry is
 * {name, role, reason}.
 */
function buildSequenceDictionary(report, options) {
  const primary = options.primary;
  const alternateColumns = (options.alternates || []).filter(a => a !== primary);
  const roles = options.roles || [];
  const roleIndex = new Map(roles.map((r, i) => [r, i]));
  const assemblyId = (options.assemblyId || "").trim();
  const warnings = [];
  const skipped = [];

  let asValue = "";
  if (assemblyId !== "") {
    if (TagValueRegex.test(assemblyId)) asValue = assemblyId;
    else warnings.push(`The assembly identifier (AS) "${assemblyId}" contains characters other than printable ASCII and was omitted.`);
  }

  // Select and order the records by role
  const candidates = [];
  report.records.forEach((record, index) => {
    const role = record[SequenceRoleColumn.key];
    if (roleIndex.has(role)) candidates.push({ record, index });
    else skipped.push({ name: recordLabel(record, primary), role, reason: `sequence role "${role}" not selected` });
  });
  if (options.sortByRole) {
    candidates.sort((a, b) =>
      roleIndex.get(a.record[SequenceRoleColumn.key]) - roleIndex.get(b.record[SequenceRoleColumn.key]) || a.index - b.index);
  }

  // Primary names, lengths, and candidate alternate names
  const sequences = [];
  const seenPrimary = new Set();
  for (const { record } of candidates) {
    const role = record[SequenceRoleColumn.key];
    const name = nameValue(record, primary);
    if (name === null) {
      const raw = record[primary];
      const reason = primary === AssignedMoleculeKey && role !== AssembledMoleculeRole && raw && raw !== MissingColumnValue
        ? `${AssignedMoleculeKey} names the parent molecule, not the sequence, for "${role}" sequences`
        : `missing a value for the primary column "${primary}"`;
      skipped.push({ name: recordLabel(record, primary), role, reason });
      continue;
    }
    if (!SequenceNameRegex.test(name)) {
      skipped.push({ name, role, reason: "not a valid SAM reference sequence name" });
      continue;
    }
    if (seenPrimary.has(name)) {
      skipped.push({ name, role, reason: "duplicate primary name" });
      continue;
    }
    const lengthText = record[SequenceLengthKey];
    const length = /^\d+$/.test(lengthText) ? Number(lengthText) : NaN;
    if (!(length >= 1 && length <= MaxSequenceLength)) {
      skipped.push({ name, role, reason: `invalid sequence length "${lengthText}"` });
      continue;
    }
    const candidateAlternates = [];
    for (const column of alternateColumns) {
      const value = nameValue(record, column);
      if (value !== null && value !== name) candidateAlternates.push([column, value]);
    }
    seenPrimary.add(name);
    sequences.push({ name, length, role, record, candidateAlternates, altByColumn: new Map(), alternates: [] });
  }

  // Resolve alternate names against the sequences that are output.  Skipping a sequence for having no alternate
  // names frees its primary name and its candidate alternates, which may give another sequence an alternate.  So
  // each round skips only the sequences with no alternates that cannot gain one by another skip, then resolves
  // again.  If every sequence without alternates is waiting on another one, the waits form cycles; then one sequence
  // is skipped from a cycle that waits on nothing outside itself (the first in output order), and the rest are
  // resolved again, since skipping it may free a name for the others.  Each round shrinks the kept set, so this
  // terminates.
  let kept = sequences;
  let altWarnings = [];
  for (;;) {
    altWarnings = resolveAlternates(kept);
    if (!options.skipMissingAlternates) break;
    const empty = kept.filter(seq => seq.alternates.length === 0);
    if (empty.length === 0) break;
    const emptySet = new Set(empty);
    const byName = new Map(kept.map(seq => [seq.name, seq]));
    // The sequences without alternates that a sequence waits on: for each valid candidate, the sequence whose primary
    // name it is and the other sequences claiming it, provided all of them also have no alternates.  A candidate
    // blocked by any sequence that has alternates stays blocked, so it adds nothing.
    const waitsOn = seq => {
      const blockers = new Set();
      for (const [, value] of seq.candidateAlternates) {
        if (!AlternateNameRegex.test(value)) continue;
        const owner = byName.get(value);
        const others = kept.filter(o => o !== seq && (o === owner || o.candidateAlternates.some(([, v]) => v === value)));
        if (others.length > 0 && others.every(o => emptySet.has(o))) others.forEach(o => blockers.add(o));
      }
      return blockers;
    };
    const waits = new Map(empty.map(seq => [seq, waitsOn(seq)]));
    let remove = empty.filter(seq => waits.get(seq).size === 0);
    if (remove.length === 0) {
      // Every sequence is waiting, so the waits form cycles.  Skip the first sequence on a cycle that waits on
      // nothing outside itself: one whose every reachable sequence can reach it back.
      const reachable = new Map();
      const reach = seq => {
        if (!reachable.has(seq)) {
          const seen = new Set();
          const stack = [...waits.get(seq)];
          while (stack.length > 0) {
            const next = stack.pop();
            if (seen.has(next)) continue;
            seen.add(next);
            for (const blocker of waits.get(next)) stack.push(blocker);
          }
          reachable.set(seq, seen);
        }
        return reachable.get(seq);
      };
      remove = [empty.find(seq => [...reach(seq)].every(other => reach(other).has(seq)))];
    }
    const removed = new Set(remove);
    for (const seq of remove) skipped.push({ name: seq.name, role: seq.role, reason: "no alternate names" });
    kept = kept.filter(seq => !removed.has(seq));
  }
  for (const warning of altWarnings) warnings.push(warning);

  // Render the dictionary
  const lines = ["@HD\tVN:1.6"];
  let totalLength = 0;
  for (const seq of kept) {
    const fields = ["@SQ", `SN:${seq.name}`, `LN:${seq.length}`];
    if (asValue !== "") fields.push(`AS:${asValue}`);
    if (seq.alternates.length > 0) fields.push(`AN:${seq.alternates.join(",")}`);
    if (options.includeReportTags) {
      const tags = [[SequenceRoleColumn.tag, seq.role]].concat(NameColumns.map(c => [c.tag, seq.record[c.key]]));
      for (const [tag, value] of tags) {
        if (value === undefined) continue; // column not in the report
        if (TagValueRegex.test(value)) fields.push(`${tag}:${value}`);
        else warnings.push(`${seq.name}: omitted tag ${tag} because its value "${value}" is empty or not printable ASCII.`);
      }
    }
    lines.push(fields.join("\t"));
    totalLength += seq.length;
  }

  return { sequences: kept, skipped, warnings, text: lines.join("\n") + "\n", totalLength, alternateColumns };
}

/**
 * Sets `alternates` and `altByColumn` on each sequence from its candidate alternates, dropping names that are not
 * valid, that are another sequence's primary name, or that are shared by more than one sequence.  Returns warnings.
 */
function resolveAlternates(sequences) {
  const warnings = [];
  const primaryNames = new Set(sequences.map(s => s.name));
  const valid = new Map(); // sequence -> [column, value][]
  const claimants = new Map(); // value -> Set of sequence names
  for (const seq of sequences) {
    const keep = [];
    for (const [column, value] of seq.candidateAlternates) {
      if (!AlternateNameRegex.test(value)) {
        warnings.push(`${seq.name}: dropped alternate name "${value}" from ${column} because it is not a valid SAM name.`);
      } else if (primaryNames.has(value)) {
        warnings.push(`${seq.name}: dropped alternate name "${value}" from ${column} because it is the primary name of another sequence.`);
      } else {
        keep.push([column, value]);
        if (!claimants.has(value)) claimants.set(value, new Set());
        claimants.get(value).add(seq.name);
      }
    }
    valid.set(seq, keep);
  }
  for (const [value, names] of claimants) {
    if (names.size > 1) {
      warnings.push(`Dropped alternate name "${value}" from ${[...names].join(", ")} because it is shared by more than one sequence.`);
    }
  }
  for (const seq of sequences) {
    seq.altByColumn = new Map();
    seq.alternates = [];
    for (const [column, value] of valid.get(seq)) {
      if (claimants.get(value).size > 1) continue;
      seq.altByColumn.set(column, value);
      if (!seq.alternates.includes(value)) seq.alternates.push(value);
    }
  }
  return warnings;
}

/** Parses a GenBank/RefSeq assembly accession, returning {prefix, digits, version} or null. */
function parseAccession(text) {
  const m = String(text).trim().match(/^(GC[AF])_(\d{9})(?:\.(\d+))?$/i);
  if (!m) return null;
  return { prefix: m[1].toUpperCase(), digits: m[2], version: m[3] === undefined ? null : Number(m[3]) };
}

/** Returns the NCBI FTP directory URL that holds all versions of the given parsed accession. */
function accessionDirectoryUrl(acc) {
  const d = acc.digits;
  return `https://ftp.ncbi.nlm.nih.gov/genomes/all/${acc.prefix}/${d.slice(0, 3)}/${d.slice(3, 6)}/${d.slice(6, 9)}/`;
}

/**
 * Given the link targets (hrefs) in an NCBI FTP accession directory listing, returns the assembly directory name
 * (e.g. `GCF_000001405.40_GRCh38.p14`) for the requested version, or for the latest version if none is given.
 * Returns null if not found.
 */
function findAssemblyDirectory(hrefs, acc) {
  const re = new RegExp(`^(?:.*/)?(${acc.prefix}_${acc.digits}\\.(\\d+)_[^/]+)/?$`);
  let best = null;
  for (const href of hrefs) {
    const m = re.exec(href);
    if (m === null) continue;
    const version = Number(m[2]);
    if (acc.version !== null && version !== acc.version) continue;
    if (best === null || version > best.version) best = { name: m[1], version };
  }
  return best === null ? null : best.name;
}

/** Returns the URL of the assembly report within an NCBI assembly directory. */
function assemblyReportUrl(acc, directoryName) {
  return `${accessionDirectoryUrl(acc)}${directoryName}/${directoryName}_assembly_report.txt`;
}

/** Formats an integer with thousands separators. */
function withThousands(value) {
  return NumberFormat.format(value);
}

if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    MissingColumnValue, NameColumns, SequenceRoleColumn, AlternateColumnOrder, DefaultAlternateColumns,
    PrimaryColumnPreference, DefaultSequenceRoles, KnownSequenceRoles, SequenceNameRegex, AlternateNameRegex,
    nameColumn, parseAssemblyReport, inferAssemblyId, rolesInReport, roleCounts, nameValue, recordLabel,
    usableNameColumns, defaultColumns, buildSequenceDictionary, parseAccession, accessionDirectoryUrl,
    findAssemblyDirectory, assemblyReportUrl, withThousands
  };
}
