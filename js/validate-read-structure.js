'use strict'

/**
 * Method for validating read structures.
 */

var AnyLengthChar = '+';

function toInt(c) { return parseInt(c, 10); }
function isDigit(d) { return !isNaN(parseInt(d, 10)); }

function invalid(msg, rs, start, end) {
	var prefix = rs.substring(0, start);
	var error  = rs.substring(start, end);
	var suffix;
	if (end === rs.length) { suffix = ''; } else { suffix = rs.substring(end, rs.length); }
	return msg + ': ' + prefix + '[' + error + ']' + suffix;
}

function validateReadStructure(readStructureString) {
	if (typeof(readStructureString) === "undefined" || readStructureString === null) {
		return {"message" : 'Read structure was empty'};
	} 
	var rs = readStructureString.replace(/^\\s*|\\s*$/g,'').toUpperCase(); // trim() not available on BaseSpace
		
	if (rs.length == 0) {
		return {"message" : 'Read structure was empty'};
	}

	var i = 0;
	var segments = [];
	while (i < rs.length) {
		// Stash the beginning position of our parsing so we can highlight what we're having trouble with
		var parsePosition = i;

		// Parse out the length segment which many be 1 or more digits or the AnyLengthChar
		var c = rs.charAt(i);
		var segLength;
		if (c === AnyLengthChar) {
			i += 1;
			segLength = null;
		}
		else if (isDigit(c)) {
			segLength = 0;
			while (i < rs.length && isDigit(rs.charAt(i))) { segLength = (segLength*10) + toInt(rs.charAt(i)); i += 1; }
		}
		else {
			return {"message" : invalid('Read structure missing length information', rs, parsePosition, parsePosition+1)};
		}

		// Parse out the operator and make a segment
		if (i === rs.length) {
			return {"message" : invalid('Read structure with invalid segment', rs, parsePosition, i)};
		}
		else {
			var code = rs.charAt(i);
			i += 1;
			if (code !== 'T' && code !== 'B' && code !== 'M' && code !== 'S') {
				return {"message" : invalid('Read structure segment had unknown type', rs, parsePosition, i)};
			}
			else if (segLength !== null && segLength <= 0) {
				return {"message" : 'Read structure contained zero length segments: ' + readStructureString};
			}
			else {
				var segment = [segLength, code];
				segments.push(segment);
			}
		}
	}

	var segmentsLength = segments.length;
	if (segmentsLength == 0) {
		return {"message" : 'Read structure was empty'};
	}
	for (i = 0; i < segmentsLength-1; i++) {
		if (segments[i][0] === null) {
			return {"message" : 'Variable length (' + AnyLengthChar + ') can only be used in the last segment: ' + readStructureString};
		}
	}

	var table = "<table><tr><th>Length</th><th>Read Type</th></tr>";
	for (i = 0; i < segmentsLength; i++) {
		var segLength = segments[i][0];
		var code = segments[i][1];

		table += "<tr>";
		if (segLength === null) {
			table += "<td>Remaining</td>";
		}
		else {
			table += "<td>" + segLength.toString() + "</td>";
		}

		table += "<td>";
		switch(code) {
			case 'T':
				table += "Template";
				break;
			case 'B':
				table += "Sample Barcode";
				break;
			case 'M':
				table += "Molecular Barcode";
				break;
			case 'S':
				table += "Skipped Bases";
				break;
			default:
				table += "Bug: Unknown Read Type '" + code + "'";
		}
		table += "</td>";
		table += "</tr>";
	}

	table += "</table>";

	return {"message" : null, "table" : table};
}
