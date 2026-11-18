'use strict'

/**
 * Method for validating sample sheets.
 */

var SampleIdKey = 'Sample_ID'.toLowerCase();
var SampleNameKey = 'Sample_Name'.toLowerCase();
var LibraryIdKey = 'Library_ID'.toLowerCase();

function sampleSheetError(line, lineNumber, message) {
	var firstLine = 'Error: on line #' + lineNumber + ': ' + message;
	if (line === null) {
		return firstLine + '.';
	}
	else if (line.length > 80) {
		return firstLine + ': ' + line.substring(0, 77) + '...';
	}
	else {
		return firstLine + ': ' + line;
	}
}
      
// since basespace doesn't support indexOf on arrays WTF
function indexOf(arr, elem) {
	for (var i = 0; i < arr.length ; i++) {
		if (arr[i] === elem) return i;
	}
	return -1;
}


function validateSampleSheet(sampleSheetString, sampleBarcodeKey) {
	sampleBarcodeKey = sampleBarcodeKey.toLowerCase();

	var lineNumber = 1;
	var lines = sampleSheetString.split(/[\r\n]+/);

	// skip to '[Data]' if exists
	var i;
	for (i = 0; i < lines.length; i++, lineNumber++) {
		if (lines[i].indexOf('[Data]') !== -1) {
			i++;
			lineNumber++;
			break;
		}
	}

	// remove trailing empty lines
	var j = lines.length - 1;
	for (; 0 <= j; j--) {
		if (lines[j].length > 0) break;
	}

	// update lines
	if (i < lines.length) lines = lines.slice(i, j+1);
	else {
		lines = lines.slice(0, j+1);
		lineNumber = 1;
	}
	i = 0;

	// validate # of remaining lines
	if (lines.length === 0) return sampleSheetError(null, lineNumber, 'no header line found');
	else if (lines.length === 1) return sampleSheetError(null, lineNumber, 'no sample data found');

	// get the header
	var header = lines[i].split(',');
	for (var j = 0; j < header.length; j++) {
		header[j] = header[j].toLowerCase();
	}
	
	// validate required keys in the header
	var hasLibraryId = indexOf(header, LibraryIdKey) !== -1;
	var keysToValidate = [SampleIdKey, SampleNameKey, sampleBarcodeKey];
	for (var j = 0; j < keysToValidate.length; j++) {
		var key = keysToValidate[j];
		if (indexOf(header, key) === -1) return sampleSheetError(lines[i], lineNumber+i, 'missing key \'' + key + '\' in the header');
	}
	i++;

	// create the samples
	var samples = []
	for (; i < lines.length; i++) {
		var data = lines[i].split(',');
		if (data.length < header.length) return sampleSheetError(lines[i], lineNumber+i, 'sample line had too few values');
		else if (data.length > header.length) return sampleSheetError(lines[i], lineNumber+i, 'sample line had too many values');
		var sample = {};
		for (var j = 0; j < header.length; j++) sample[header[j]] = data[j];
		// Default Library_ID to Sample_ID for validation. NB: should have already validated that SampleIdKey exists in the header
		if (!hasLibraryId) sample[LibraryIdKey] = sample[SampleIdKey];
		samples.push(sample);
	}

	// validate that sample identifiers are unique
	var sampleIds = [];
	for (var j = 0; j < samples.length; j++) {
		var sampleId = samples[j][SampleIdKey];
		if (indexOf(sampleIds, sampleId) !== -1) return sampleSheetError(lines[j+1], lineNumber+j+1, 'key \'' + SampleIdKey + '\' was not unique');
		sampleIds.push(sampleId);
	}

	// validate that sample name and library identifiers in combination are unique
	var sampleNameAndLibraryIds = [];
	for (var j = 0; j < samples.length; j++) {
		var sampleName = samples[j][SampleNameKey];
		var libraryId = samples[j][LibraryIdKey];
		var sampleNameAndLibraryId = sampleName + '-' + libraryId;
		if (indexOf(sampleNameAndLibraryIds, sampleNameAndLibraryId) !== -1) {
			if (hasLibraryId) {
				return sampleSheetError(lines[j+1], lineNumber+j+1, 'combination of keys \'' + SampleNameKey + '\' and \'' + LibraryIdKey + '\' were not unique');
			}
			else {
				// Should never happen, since Sample_ID is unique
				return sampleSheetError(lines[j+1], lineNumber+j+1, 'combination of keys \'' + SampleNameKey + '\' and \'' + SampleId + '\' were not unique');
			}
		}
		sampleNameAndLibraryIds.push(sampleNameAndLibraryId);
	}

	// validate that the barcode column contains only DNA bases and are unique
	var sampleBarcodes = [];
	for (var j = 0; j < samples.length; j++) {
		var sampleBarcode = samples[j][sampleBarcodeKey];
		var matches = sampleBarcode.match(/^[ACGTNacgtn]+$/);
		if (matches === null) return sampleSheetError(lines[j+1], lineNumber+j+1, 'sample barcode column \'' + sampleBarcodeKey + '\' did not contain only DNA bases');
		if (indexOf(sampleBarcodes, sampleBarcode) !== -1) return sampleSheetError(lines[j+1], lineNumber+j+1, 'sample barcode \'' + sampleBarcode + '\' was not unique');
		sampleBarcodes.push(sampleBarcode);
	}

	return null;
}
