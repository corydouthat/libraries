/*
fileHelp.h
Simple file helpers
by: Connor Douthat
3/19/2015
*/
// Used by permission

#ifndef _CDSTD_FILEHELP_H_
#define _CDSTD_FILEHELP_H_

char *getFileContents(const char *filepath, unsigned int *sizeOut = NULL)
{
	// Cory Douthat - Changed to fopen_s
	// Cory Douthat - TODO: check return from fopen_s
	FILE *file;
	fopen_s(&file, filepath, "rb");
	if(!file) return NULL;

	fseek(file, 0, SEEK_END);
	unsigned int size = ftell(file);
	if(sizeOut) *sizeOut = size;
	rewind(file);

	char *contents = (char*)malloc(size + 1);
	contents[size] = 0;
	if(fread(contents, 1, size, file) != size)
	{
		free(contents);
		contents = NULL;
	}

	fclose(file);
	return contents;
}
bool setFileContents(const char *filepath, char *data, unsigned int len = 0)
{
	if(!len) len = strlen(data);

	// Cory Douthat - Changed to fopen_s
	// Cory Douthat - TODO: check return from fopen_s
	FILE *file;
	fopen_s(&file, filepath, "wb");
	if(!file) return false;

	unsigned int written = fwrite(data, 1, len, file);
	fclose(file);

	return (written == len);
}

#endif
