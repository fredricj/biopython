#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define BYTESWAP(x) { \
    uint32_t t = x; \
    x = (t >> 24) | ((t & 0x00FF0000) >> 8) | ((t & 0x0000FF00) << 8) | (t << 24); \
}

int main(int argc, char *argv[])
{
    int isByteSwapped = 0;
    int32_t signature;
    int32_t version;
    uint32_t sequenceCount;
    uint32_t reserved;
    FILE *ptr;
    uint32_t i, j;
    uint8_t nameSize;
    char name[256];
    uint32_t offset;
    uint32_t dnaSize;
    uint32_t nBlockCount;
    uint32_t* nBlockStarts;
    uint32_t* nBlockSizes;
    uint32_t maskBlockCount;
    uint32_t* maskBlockStarts;
    uint32_t* maskBlockSizes;
    unsigned char* packedDna;
    size_t n;
    if (argc != 2) {
        printf("missing input file name\n");
        return 0;
    }
    ptr = fopen(argv[1], "rt");  /* r for read, b for binary */
    fread(&signature, 4, 1, ptr);
    switch (signature) {
        case 0x1A412743: printf("signature: same architecture\n"); break;
        case 0x4327411a: printf("signature: byte-swapped\n"); isByteSwapped = 1; break;
        default:
            printf("Unknown signature %X\n", signature);
            return 0;
    }
    fread(&version, 4, 1, ptr);
    if (version != 0) {
        printf("Found non-zero file version %u; aborting\n", version);
        return 0;
    }
    fread(&sequenceCount, 4, 1, ptr);
    if (isByteSwapped) BYTESWAP(sequenceCount);
    printf("Number of sequences is %d\n", sequenceCount);
    fread(&reserved, 4, 1, ptr);
    if (reserved != 0) {
        printf("Found non-zero reserved field %u; aborting\n", reserved);
        return 0;
    }
    for (i = 0; i < sequenceCount; i++) {
        fread(&nameSize, 1, 1, ptr);
        fread(name, nameSize, 1, ptr);
        /* nameSize <= 255 */
        name[nameSize] = '\0';
        printf("sequence %u name is %s nameSize is %u\n", i, name, nameSize);
        fread(&offset, 4, 1, ptr);
        if (isByteSwapped) BYTESWAP(offset);
        printf("sequence %u offset is %u\n", i, offset);
    }
    for (i = 0; i < sequenceCount; i++) {
        fread(&dnaSize, 4, 1, ptr);
        if (isByteSwapped) BYTESWAP(dnaSize);
        printf("sequence %u dnaSize is %u\n", i, dnaSize);
        fread(&nBlockCount, 4, 1, ptr);
        if (isByteSwapped) BYTESWAP(nBlockCount);
        printf("sequence %u number of blocks of Ns is %u\n", i, nBlockCount);
        nBlockStarts = malloc(nBlockCount*sizeof(uint32_t));
        nBlockSizes = malloc(nBlockCount*sizeof(uint32_t));
        fread(nBlockStarts, 4, nBlockCount, ptr);
        fread(nBlockSizes, 4, nBlockCount, ptr);
        if (isByteSwapped) {
            for (j = 0; j < nBlockCount; j++) {
                BYTESWAP(nBlockStarts[j]);
                BYTESWAP(nBlockSizes[j]);
            }
        }
        for (j = 0; j < nBlockCount; j++) printf("nBlock start = %u size = %u\n", nBlockStarts[j], nBlockSizes[j]);
        fread(&maskBlockCount, 4, 1, ptr);
        if (isByteSwapped) BYTESWAP(maskBlockCount);
        printf("sequence %u number of masked blocks is %u\n", i, maskBlockCount);
        maskBlockStarts = malloc(maskBlockCount*sizeof(uint32_t));
        maskBlockSizes = malloc(maskBlockCount*sizeof(uint32_t));
        fread(maskBlockStarts, 4, maskBlockCount, ptr);
        fread(maskBlockSizes, 4, maskBlockCount, ptr);
        if (isByteSwapped) {
            for (j = 0; j < maskBlockCount; j++) {
                BYTESWAP(maskBlockStarts[j]);
                BYTESWAP(maskBlockSizes[j]);
            }
        }
        for (j = 0; j < maskBlockCount; j++) printf("nMask start = %u size = %u\n", maskBlockStarts[j], maskBlockSizes[j]);
        fread(&reserved, 4, 1, ptr);
        if (reserved != 0) {
            printf("Found non-zero reserved field %u in sequence record %i; aborting\n", reserved, i);
            return 0;
        }
        n = (dnaSize + 3) / 4;
        packedDna = malloc(n*sizeof(unsigned char));
        fread(packedDna, 1, n, ptr);
        printf("sequence %u, length %u, starts with %x %x %x %x, ends with %x %x %x %x\n", i, dnaSize, packedDna[0], packedDna[1], packedDna[2], packedDna[3], packedDna[n-4], packedDna[n-3], packedDna[n-2], packedDna[n-1]);
    }
    return 0;
}
