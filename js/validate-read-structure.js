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
	var rs = readStructureString.trim().toUpperCase();
		
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
		else if (isDigit(c) || c === '-') {
			var sign = 1;
			if (c === '-') {
				var sign = -1;
				i++;
			}

			segLength = 0;
			while (i < rs.length && isDigit(rs.charAt(i))) { segLength = (segLength*10) + toInt(rs.charAt(i)); i += 1; }
			segLength *= sign;
			
			if (segLength < 0) {
				return {"message" : invalid('Read structure contained a negative length segment', rs, parsePosition, i)};
			}
			else if (segLength == 0) {
				return {"message" : invalid('Read structure contained a zero length segment', rs, parsePosition, i)};
			}
		}
		else {
			return {"message" : invalid('Read structure segment missing length information', rs, parsePosition, parsePosition+1)};
		}

		// Parse out the operator and make a segment
		if (i === rs.length) {
			return {"message" : invalid('Read structure had an invalid segment', rs, parsePosition, i)};
		}
		else {
			var code = rs.charAt(i);
			i += 1;
			if (code !== 'T' && code !== 'B' && code !== 'M' && code !== 'S') {
				return {"message" : invalid('Read structure segment had unknown type', rs, parsePosition, i)};
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

	return {"message" : null, "segments" : segments};
}

function validateReadStructures(readStructuresString) {
	var readStructures = readStructuresString.trim().split(/\s+/);
	var i = 0;
	var table = "";
	if (readStructures.length == 1) {
		table = "<table><tr><th>Length</th><th>Read Type</th></tr>";
	}
	else {
		table = "<table><tr><th>Read Structure</th><th>Length</th><th>Read Type</th></tr>";
	}
	for (i = 0; i < readStructures.length; i++) {
		var result = validateReadStructure(readStructures[i]);
		if (result["message"] !== null) {
			return result;
		}

		var segments = result["segments"];
		var segmentsLength = segments.length;
		var j = 0;
		for (j = 0; j < segmentsLength; j++) {
			var segLength = segments[j][0];
			var code = segments[j][1];

			table += "<tr>";
			if (readStructures.length > 1) {
				table += "<td>" + (i+1) + "</td>";
			}
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
	}

	table += "</table>";

	return {"message": null, "table" : table};
}
