#include "Python.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>


#define BYTESWAP(x) { \
    uint32_t t = x; \
    x = (t >> 24) | ((t & 0x00FF0000) >> 8) | ((t & 0x0000FF00) << 8) | (t << 24); \
}

static const char bases[][4] = {"TTTT",  /* 00 00 00 00 */
                                "TTTC",  /* 00 00 00 01 */
                                "TTTA",  /* 00 00 00 10 */
                                "TTTG",  /* 00 00 00 11 */
                                "TTCT",  /* 00 00 01 00 */
                                "TTCC",  /* 00 00 01 01 */
                                "TTCA",  /* 00 00 01 10 */
                                "TTCG",  /* 00 00 01 11 */
                                "TTAT",  /* 00 00 10 00 */
                                "TTAC",  /* 00 00 10 01 */
                                "TTAA",  /* 00 00 10 10 */
                                "TTAG",  /* 00 00 10 11 */
                                "TTGT",  /* 00 00 11 00 */
                                "TTGC",  /* 00 00 11 01 */
                                "TTGA",  /* 00 00 11 10 */
                                "TTGG",  /* 00 00 11 11 */
                                "TCTT",  /* 00 01 00 00 */
                                "TCTC",  /* 00 01 00 01 */
                                "TCTA",  /* 00 01 00 10 */
                                "TCTG",  /* 00 01 00 11 */
                                "TCCT",  /* 00 01 01 00 */
                                "TCCC",  /* 00 01 01 01 */
                                "TCCA",  /* 00 01 01 10 */
                                "TCCG",  /* 00 01 01 11 */
                                "TCAT",  /* 00 01 10 00 */
                                "TCAC",  /* 00 01 10 01 */
                                "TCAA",  /* 00 01 10 10 */
                                "TCAG",  /* 00 01 10 11 */
                                "TCGT",  /* 00 01 11 00 */
                                "TCGC",  /* 00 01 11 01 */
                                "TCGA",  /* 00 01 11 10 */
                                "TCGG",  /* 00 01 11 11 */
                                "TATT",  /* 00 10 00 00 */
                                "TATC",  /* 00 10 00 01 */
                                "TATA",  /* 00 10 00 10 */
                                "TATG",  /* 00 10 00 11 */
                                "TACT",  /* 00 10 01 00 */
                                "TACC",  /* 00 10 01 01 */
                                "TACA",  /* 00 10 01 10 */
                                "TACG",  /* 00 10 01 11 */
                                "TAAT",  /* 00 10 10 00 */
                                "TAAC",  /* 00 10 10 01 */
                                "TAAA",  /* 00 10 10 10 */
                                "TAAG",  /* 00 10 10 11 */
                                "TAGT",  /* 00 10 11 00 */
                                "TAGC",  /* 00 10 11 01 */
                                "TAGA",  /* 00 10 11 10 */
                                "TAGG",  /* 00 10 11 11 */
                                "TGTT",  /* 00 11 00 00 */
                                "TGTC",  /* 00 11 00 01 */
                                "TGTA",  /* 00 11 00 10 */
                                "TGTG",  /* 00 11 00 11 */
                                "TGCT",  /* 00 11 01 00 */
                                "TGCC",  /* 00 11 01 01 */
                                "TGCA",  /* 00 11 01 10 */
                                "TGCG",  /* 00 11 01 11 */
                                "TGAT",  /* 00 11 10 00 */
                                "TGAC",  /* 00 11 10 01 */
                                "TGAA",  /* 00 11 10 10 */
                                "TGAG",  /* 00 11 10 11 */
                                "TGGT",  /* 00 11 11 00 */
                                "TGGC",  /* 00 11 11 01 */
                                "TGGA",  /* 00 11 11 10 */
                                "TGGG",  /* 00 11 11 11 */
                                "CTTT",  /* 01 00 00 00 */
                                "CTTC",  /* 01 00 00 01 */
                                "CTTA",  /* 01 00 00 10 */
                                "CTTG",  /* 01 00 00 11 */
                                "CTCT",  /* 01 00 01 00 */
                                "CTCC",  /* 01 00 01 01 */
                                "CTCA",  /* 01 00 01 10 */
                                "CTCG",  /* 01 00 01 11 */
                                "CTAT",  /* 01 00 10 00 */
                                "CTAC",  /* 01 00 10 01 */
                                "CTAA",  /* 01 00 10 10 */
                                "CTAG",  /* 01 00 10 11 */
                                "CTGT",  /* 01 00 11 00 */
                                "CTGC",  /* 01 00 11 01 */
                                "CTGA",  /* 01 00 11 10 */
                                "CTGG",  /* 01 00 11 11 */
                                "CCTT",  /* 01 01 00 00 */
                                "CCTC",  /* 01 01 00 01 */
                                "CCTA",  /* 01 01 00 10 */
                                "CCTG",  /* 01 01 00 11 */
                                "CCCT",  /* 01 01 01 00 */
                                "CCCC",  /* 01 01 01 01 */
                                "CCCA",  /* 01 01 01 10 */
                                "CCCG",  /* 01 01 01 11 */
                                "CCAT",  /* 01 01 10 00 */
                                "CCAC",  /* 01 01 10 01 */
                                "CCAA",  /* 01 01 10 10 */
                                "CCAG",  /* 01 01 10 11 */
                                "CCGT",  /* 01 01 11 00 */
                                "CCGC",  /* 01 01 11 01 */
                                "CCGA",  /* 01 01 11 10 */
                                "CCGG",  /* 01 01 11 11 */
                                "CATT",  /* 01 10 00 00 */
                                "CATC",  /* 01 10 00 01 */
                                "CATA",  /* 01 10 00 10 */
                                "CATG",  /* 01 10 00 11 */
                                "CACT",  /* 01 10 01 00 */
                                "CACC",  /* 01 10 01 01 */
                                "CACA",  /* 01 10 01 10 */
                                "CACG",  /* 01 10 01 11 */
                                "CAAT",  /* 01 10 10 00 */
                                "CAAC",  /* 01 10 10 01 */
                                "CAAA",  /* 01 10 10 10 */
                                "CAAG",  /* 01 10 10 11 */
                                "CAGT",  /* 01 10 11 00 */
                                "CAGC",  /* 01 10 11 01 */
                                "CAGA",  /* 01 10 11 10 */
                                "CAGG",  /* 01 10 11 11 */
                                "CGTT",  /* 01 11 00 00 */
                                "CGTC",  /* 01 11 00 01 */
                                "CGTA",  /* 01 11 00 10 */
                                "CGTG",  /* 01 11 00 11 */
                                "CGCT",  /* 01 11 01 00 */
                                "CGCC",  /* 01 11 01 01 */
                                "CGCA",  /* 01 11 01 10 */
                                "CGCG",  /* 01 11 01 11 */
                                "CGAT",  /* 01 11 10 00 */
                                "CGAC",  /* 01 11 10 01 */
                                "CGAA",  /* 01 11 10 10 */
                                "CGAG",  /* 01 11 10 11 */
                                "CGGT",  /* 01 11 11 00 */
                                "CGGC",  /* 01 11 11 01 */
                                "CGGA",  /* 01 11 11 10 */
                                "CGGG",  /* 01 11 11 11 */
                                "ATTT",  /* 10 00 00 00 */
                                "ATTC",  /* 10 00 00 01 */
                                "ATTA",  /* 10 00 00 10 */
                                "ATTG",  /* 10 00 00 11 */
                                "ATCT",  /* 10 00 01 00 */
                                "ATCC",  /* 10 00 01 01 */
                                "ATCA",  /* 10 00 01 10 */
                                "ATCG",  /* 10 00 01 11 */
                                "ATAT",  /* 10 00 10 00 */
                                "ATAC",  /* 10 00 10 01 */
                                "ATAA",  /* 10 00 10 10 */
                                "ATAG",  /* 10 00 10 11 */
                                "ATGT",  /* 10 00 11 00 */
                                "ATGC",  /* 10 00 11 01 */
                                "ATGA",  /* 10 00 11 10 */
                                "ATGG",  /* 10 00 11 11 */
                                "ACTT",  /* 10 01 00 00 */
                                "ACTC",  /* 10 01 00 01 */
                                "ACTA",  /* 10 01 00 10 */
                                "ACTG",  /* 10 01 00 11 */
                                "ACCT",  /* 10 01 01 00 */
                                "ACCC",  /* 10 01 01 01 */
                                "ACCA",  /* 10 01 01 10 */
                                "ACCG",  /* 10 01 01 11 */
                                "ACAT",  /* 10 01 10 00 */
                                "ACAC",  /* 10 01 10 01 */
                                "ACAA",  /* 10 01 10 10 */
                                "ACAG",  /* 10 01 10 11 */
                                "ACGT",  /* 10 01 11 00 */
                                "ACGC",  /* 10 01 11 01 */
                                "ACGA",  /* 10 01 11 10 */
                                "ACGG",  /* 10 01 11 11 */
                                "AATT",  /* 10 10 00 00 */
                                "AATC",  /* 10 10 00 01 */
                                "AATA",  /* 10 10 00 10 */
                                "AATG",  /* 10 10 00 11 */
                                "AACT",  /* 10 10 01 00 */
                                "AACC",  /* 10 10 01 01 */
                                "AACA",  /* 10 10 01 10 */
                                "AACG",  /* 10 10 01 11 */
                                "AAAT",  /* 10 10 10 00 */
                                "AAAC",  /* 10 10 10 01 */
                                "AAAA",  /* 10 10 10 10 */
                                "AAAG",  /* 10 10 10 11 */
                                "AAGT",  /* 10 10 11 00 */
                                "AAGC",  /* 10 10 11 01 */
                                "AAGA",  /* 10 10 11 10 */
                                "AAGG",  /* 10 10 11 11 */
                                "AGTT",  /* 10 11 00 00 */
                                "AGTC",  /* 10 11 00 01 */
                                "AGTA",  /* 10 11 00 10 */
                                "AGTG",  /* 10 11 00 11 */
                                "AGCT",  /* 10 11 01 00 */
                                "AGCC",  /* 10 11 01 01 */
                                "AGCA",  /* 10 11 01 10 */
                                "AGCG",  /* 10 11 01 11 */
                                "AGAT",  /* 10 11 10 00 */
                                "AGAC",  /* 10 11 10 01 */
                                "AGAA",  /* 10 11 10 10 */
                                "AGAG",  /* 10 11 10 11 */
                                "AGGT",  /* 10 11 11 00 */
                                "AGGC",  /* 10 11 11 01 */
                                "AGGA",  /* 10 11 11 10 */
                                "AGGG",  /* 10 11 11 11 */
                                "GTTT",  /* 11 00 00 00 */
                                "GTTC",  /* 11 00 00 01 */
                                "GTTA",  /* 11 00 00 10 */
                                "GTTG",  /* 11 00 00 11 */
                                "GTCT",  /* 11 00 01 00 */
                                "GTCC",  /* 11 00 01 01 */
                                "GTCA",  /* 11 00 01 10 */
                                "GTCG",  /* 11 00 01 11 */
                                "GTAT",  /* 11 00 10 00 */
                                "GTAC",  /* 11 00 10 01 */
                                "GTAA",  /* 11 00 10 10 */
                                "GTAG",  /* 11 00 10 11 */
                                "GTGT",  /* 11 00 11 00 */
                                "GTGC",  /* 11 00 11 01 */
                                "GTGA",  /* 11 00 11 10 */
                                "GTGG",  /* 11 00 11 11 */
                                "GCTT",  /* 11 01 00 00 */
                                "GCTC",  /* 11 01 00 01 */
                                "GCTA",  /* 11 01 00 10 */
                                "GCTG",  /* 11 01 00 11 */
                                "GCCT",  /* 11 01 01 00 */
                                "GCCC",  /* 11 01 01 01 */
                                "GCCA",  /* 11 01 01 10 */
                                "GCCG",  /* 11 01 01 11 */
                                "GCAT",  /* 11 01 10 00 */
                                "GCAC",  /* 11 01 10 01 */
                                "GCAA",  /* 11 01 10 10 */
                                "GCAG",  /* 11 01 10 11 */
                                "GCGT",  /* 11 01 11 00 */
                                "GCGC",  /* 11 01 11 01 */
                                "GCGA",  /* 11 01 11 10 */
                                "GCGG",  /* 11 01 11 11 */
                                "GATT",  /* 11 10 00 00 */
                                "GATC",  /* 11 10 00 01 */
                                "GATA",  /* 11 10 00 10 */
                                "GATG",  /* 11 10 00 11 */
                                "GACT",  /* 11 10 01 00 */
                                "GACC",  /* 11 10 01 01 */
                                "GACA",  /* 11 10 01 10 */
                                "GACG",  /* 11 10 01 11 */
                                "GAAT",  /* 11 10 10 00 */
                                "GAAC",  /* 11 10 10 01 */
                                "GAAA",  /* 11 10 10 10 */
                                "GAAG",  /* 11 10 10 11 */
                                "GAGT",  /* 11 10 11 00 */
                                "GAGC",  /* 11 10 11 01 */
                                "GAGA",  /* 11 10 11 10 */
                                "GAGG",  /* 11 10 11 11 */
                                "GGTT",  /* 11 11 00 00 */
                                "GGTC",  /* 11 11 00 01 */
                                "GGTA",  /* 11 11 00 10 */
                                "GGTG",  /* 11 11 00 11 */
                                "GGCT",  /* 11 11 01 00 */
                                "GGCC",  /* 11 11 01 01 */
                                "GGCA",  /* 11 11 01 10 */
                                "GGCG",  /* 11 11 01 11 */
                                "GGAT",  /* 11 11 10 00 */
                                "GGAC",  /* 11 11 10 01 */
                                "GGAA",  /* 11 11 10 10 */
                                "GGAG",  /* 11 11 10 11 */
                                "GGGT",  /* 11 11 11 00 */
                                "GGGC",  /* 11 11 11 01 */
                                "GGGA",  /* 11 11 11 10 */
                                "GGGG",  /* 11 11 11 11 */
                               };

static int
safe_read(int fd, ssize_t size, void* pointer, const char variable[])
{
    const ssize_t n = read(fd, pointer, size);
    if (n == size) return 1;
    switch (n) {
        case -1:
            PyErr_Format(PyExc_RuntimeError,
                "error occurred trying to read %s (errno = %d)",
                variable, errno);
            break;
        case 0:
        default:
            PyErr_Format(PyExc_RuntimeError,
                "unexpected end of file while reading %s", variable);
            break;
    }
    return 0;
}

static char* unpack(int fd, uint32_t start, uint32_t end) {
    uint32_t i;
    const uint32_t size = end - start;
    const uint32_t byteStart = start / 4;
    const uint32_t byteEnd = (end + 3) / 4;
    const uint32_t byteSize = byteEnd - byteStart;
    unsigned char* bytes = malloc(byteSize*sizeof(unsigned char));
    char* sequence = malloc((size+1)*sizeof(char));
    if (lseek(fd, byteStart, SEEK_CUR) == -1) {
        PyErr_SetString(PyExc_RuntimeError, "failed to seek in sequence");
        return NULL;
    }
    if (!safe_read(fd, byteSize, bytes, "sequence data")) return NULL;
    start -= byteStart * 4;
    if (byteStart + 1 == byteEnd) {
        /* one byte only */
        memcpy(sequence, &(bases[*bytes][start]), size);
    }
    else {
        end -= byteEnd * 4;
        /* end is now a negative number equal to the distance to the byte end */
        memcpy(sequence, &(bases[*bytes][start]), 4 - start);
        bytes++;
        sequence += (4 - start);
        for (i = byteStart+1; i < byteEnd-1; i++, bytes++, sequence += 4)
            memcpy(sequence, bases[*bytes], 4);
        memcpy(sequence, bases[*bytes], end + 4);
        bytes++;
        bytes -= byteSize;
        sequence += (end + 4);
        sequence -= size;
    }
    sequence[size] = '\0';
    free(bytes);
    return sequence;
}

static void
applyNs(char sequence[], uint32_t start, uint32_t end,
        uint32_t nBlockCount, uint32_t nBlockStarts[], uint32_t nBlockSizes[])
{
    uint32_t i;
    for (i = 0; i < nBlockCount; i++) {
        uint32_t nBlockStart = nBlockStarts[i];
        uint32_t nBlockEnd = nBlockStart + nBlockSizes[i];
        if (nBlockEnd < start) continue;
        if (end < nBlockStart) break;
        if (nBlockStart < start) nBlockStart = start;
        if (end < nBlockEnd) nBlockEnd = end;
        memset(sequence + nBlockStart - start, 'N', nBlockEnd - nBlockStart);
    }
}

static void
applyMask(char sequence[], uint32_t start, uint32_t end,
        uint32_t maskBlockCount, uint32_t maskBlockStarts[], uint32_t maskBlockSizes[])
{
    uint32_t i, j;
    for (i = 0; i < maskBlockCount; i++) {
        uint32_t maskBlockStart = maskBlockStarts[i];
        uint32_t maskBlockEnd = maskBlockStart + maskBlockSizes[i];
        if (maskBlockEnd < start) continue;
        if (end < maskBlockStart) break;
        if (maskBlockStart < start) maskBlockStart = start;
        if (end < maskBlockEnd) maskBlockEnd = end;
        for (j = maskBlockStart - start; j < maskBlockEnd - start; j++)
            sequence[j] += 32;
    }
}

static char TwoBitIterator__doc__[] = "Lazy parser of a TwoBit file";

typedef struct {
    PyObject_HEAD
    PyObject* fileobj;
    int fd;
    uint32_t offset;
    uint32_t dnaSize;
} TwoBitSequence;

static void
TwoBitSequence_dealloc(TwoBitSequence* self)
{
    Py_DECREF(self->fileobj);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject*
TwoBitSequence_str(TwoBitSequence* self)
{
    char buffer[128];
    sprintf(buffer, "TwoBit sequence of length %u", self->dnaSize);
    return PyUnicode_FromString(buffer);
}

static int
TwoBitSequence_length(TwoBitSequence *self)
{
    return self->dnaSize;
}

static PyObject*
TwoBitSequence_subscript(TwoBitSequence* self, PyObject* item)
{
    char buffer[256];
    const uint32_t n = self->dnaSize;
    if (PyIndex_Check(item)) {
        Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);
        if (i == -1 && PyErr_Occurred())
            return NULL;
        if (i < 0)
            i += n;
        if (i < 0 || i >= n) {
            PyErr_SetString(PyExc_IndexError,
                            "TwoBitSequence index out of range");
            return NULL;
        }
        sprintf(buffer, "item %lu", i);
        return PyUnicode_FromString(buffer);
    }
    else if (PySlice_Check(item)) {
        Py_ssize_t i, j;
        Py_ssize_t start, stop, step, slicelength;
        if (PySlice_GetIndicesEx(item, n, &start, &stop, &step,
                                 &slicelength) == -1) return NULL;
        if (slicelength == 0) return PyList_New(0);
        else {
            PyObject* text;
            PyObject* result = PyList_New(slicelength);
            if (!result) return PyErr_NoMemory();
            for (i = 0, j = start; i < slicelength; i++, j += step) {
                sprintf(buffer, "item %lu", j);
                text = PyUnicode_FromString(buffer);
                PyList_SET_ITEM(result, i, text);
            }
            return result;
        }
    }
    else {
        PyErr_Format(PyExc_TypeError,
                     "TwoBitSequence indices must be integers, not %.200s",
                     Py_TYPE(item)->tp_name);
        return NULL;
    }
}

static PyMappingMethods TwoBitSequence_mapping = {
    (lenfunc)TwoBitSequence_length,           /* mp_length */
    (binaryfunc)TwoBitSequence_subscript,     /* mp_subscript */
};

static char TwoBitSequence_doc[] =
"TwoBitSequence objects store the information in TwoBit files needed to read\n"
"their sequence contents.";

static PyTypeObject TwoBitSequenceType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_TwoBitIO.TwoBitSequence",                  /* tp_name */
    sizeof(TwoBitSequence),                      /* tp_basicsize */
    0,                                           /* tp_itemsize */
    (destructor)TwoBitSequence_dealloc,          /* tp_dealloc */
    0,                                           /* tp_print */
    0,                                           /* tp_getattr */
    0,                                           /* tp_setattr */
    0,                                           /* tp_compare */
    0,                                           /* tp_repr */
    0,                                           /* tp_as_number */
    0,                                           /* tp_as_sequence */
    &TwoBitSequence_mapping,                     /* tp_as_mapping */
    0,                                           /* tp_hash */
    0,                                           /* tp_call */
    (reprfunc)TwoBitSequence_str,                /* tp_str */
    0,                                           /* tp_getattro */
    0,                                           /* tp_setattro */
    0,                                           /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,    /* tp_flags*/
    TwoBitSequence_doc,                          /* tp_doc */
};

static PyObject*
TwoBitIterator(PyObject* self, PyObject* args, PyObject* keywords)
{
    uint32_t start;
    uint32_t end;
    int isByteSwapped = 0;
    uint32_t signature;
    uint32_t version;
    uint32_t sequenceCount;
    uint32_t reserved;
    uint32_t i, j;
    uint8_t nameSize;
    char name[256];
    uint32_t offset;
    uint32_t* offsets;
    uint32_t dnaSize;
    uint32_t nBlockCount;
    uint32_t* nBlockStarts;
    uint32_t* nBlockSizes;
    uint32_t maskBlockCount;
    uint32_t* maskBlockStarts;
    uint32_t* maskBlockSizes;
    char* sequence;
    PyObject* s;
    PyObject* sequence_tuple;
    PyObject* tuple;
    PyObject* fileobj;
    int fd;
    TwoBitSequence* tbs;

    static char* kwlist[] = {"handle",
                             "start",
                             "end",
                              NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywords, "OII", kwlist,
                                     &fileobj, &start, &end)) return NULL;

    fd = PyObject_AsFileDescriptor(fileobj);
    if (fd == -1) return NULL;

    if (!safe_read(fd, sizeof(uint32_t), &signature, "signature")) return NULL;
    switch (signature) {
        case 0x1A412743: break;
        case 0x4327411a: isByteSwapped = 1; break;
        default:
            PyErr_Format(PyExc_RuntimeError,
                         "Unknown signature %X", signature);
            return 0;
    }
    if (!safe_read(fd, sizeof(uint32_t), &version, "version")) return NULL;
    if (version != 0) {
        PyErr_Format(PyExc_RuntimeError,
                     "Found non-zero file version %u; aborting", version);
        return 0;
    }
    if (!safe_read(fd, sizeof(uint32_t), &sequenceCount, "sequenceCount")) return NULL;
    if (isByteSwapped) BYTESWAP(sequenceCount);
    if (!safe_read(fd, sizeof(uint32_t), &reserved, "reserved field")) return NULL;
    if (reserved != 0) {
        PyErr_SetString(PyExc_RuntimeError,
                        "Found non-zero reserved field; aborting");
        return 0;
    }
    offsets = malloc(sequenceCount*sizeof(uint32_t));
    for (i = 0; i < sequenceCount; i++) {
        if (!safe_read(fd, sizeof(uint8_t), &nameSize, "nameSize")) return NULL;
        if (!safe_read(fd, nameSize, name, "name")) return NULL;
        if (!safe_read(fd, sizeof(uint32_t), &offset, "offset")) return NULL;
        if (isByteSwapped) BYTESWAP(offset);
        name[nameSize] = '\0';  /* nameSize <= 255 */
        offsets[i] = offset;
    }
    sequence_tuple = PyTuple_New(sequenceCount);
    tuple = PyTuple_New(sequenceCount);
    for (i = 0; i < sequenceCount; i++) {
        if (lseek(fd, offsets[i], SEEK_SET) == -1) {
            PyErr_SetString(PyExc_RuntimeError, "failed to seek in file");
            return NULL;
        }
        if (!safe_read(fd, sizeof(uint32_t), &dnaSize, "dnaSize")) return NULL;
        if (isByteSwapped) BYTESWAP(dnaSize);
        if (!safe_read(fd, sizeof(uint32_t), &nBlockCount, "nBlockCount")) return NULL;
        if (isByteSwapped) BYTESWAP(nBlockCount);
        nBlockStarts = malloc(nBlockCount*sizeof(uint32_t));
        nBlockSizes = malloc(nBlockCount*sizeof(uint32_t));
        for (j = 0; j < nBlockCount; j++) {
            if (!safe_read(fd, sizeof(uint32_t),
                           nBlockStarts+j, "nBlockStarts")) return NULL;
            if (isByteSwapped) BYTESWAP(nBlockStarts[j]);
        }
        for (j = 0; j < nBlockCount; j++) {
            if (!safe_read(fd, sizeof(uint32_t),
                           nBlockSizes+j, "nBlockSizes")) return NULL;
            if (isByteSwapped) BYTESWAP(nBlockSizes[j]);
        }
        if (!safe_read(fd, sizeof(uint32_t), &maskBlockCount, "maskBlockCount")) return NULL;
        if (isByteSwapped) BYTESWAP(maskBlockCount);
        maskBlockStarts = malloc(maskBlockCount*sizeof(uint32_t));
        maskBlockSizes = malloc(maskBlockCount*sizeof(uint32_t));
        for (j = 0; j < maskBlockCount; j++) {
            if (!safe_read(fd, sizeof(uint32_t),
                           maskBlockStarts+j, "maskBlockStarts")) return NULL;
            if (isByteSwapped) BYTESWAP(maskBlockStarts[j]);
        }
        for (j = 0; j < maskBlockCount; j++) {
            if (!safe_read(fd, sizeof(uint32_t),
                           maskBlockSizes+j, "maskBlockSizes")) return NULL;
            if (isByteSwapped) BYTESWAP(maskBlockSizes[j]);
        }
        if (!safe_read(fd, sizeof(uint32_t),
                       &reserved, "reserved")) return NULL;
        if (reserved != 0) {
            PyErr_Format(PyExc_RuntimeError,
                         "Found non-zero reserved field %u in sequence "
                         "record %i; aborting\n", reserved, i);
            free(nBlockStarts);
            free(nBlockSizes);
            free(maskBlockStarts);
            free(maskBlockSizes);
            free(offsets);
            return 0;
        }
        sequence = unpack(fd, start, end);
        applyNs(sequence, start, end, nBlockCount, nBlockStarts, nBlockSizes);
        applyMask(sequence, start, end,
                  maskBlockCount, maskBlockStarts, maskBlockSizes);
        s = PyUnicode_FromString(sequence);
        PyTuple_SET_ITEM(sequence_tuple, i, s);

        tbs = (TwoBitSequence*)PyType_GenericAlloc(&TwoBitSequenceType, 0);
        if (!tbs) return NULL;
        Py_INCREF(fileobj);
        tbs->fileobj = fileobj;
        tbs->fd = fd;
        tbs->offset = offset;
        tbs->dnaSize = dnaSize;
        PyTuple_SET_ITEM(tuple, i, (PyObject*)tbs);

        free(sequence);
        free(nBlockStarts);
        free(nBlockSizes);
        free(maskBlockStarts);
        free(maskBlockSizes);
    }
    free(offsets);
    return Py_BuildValue("OOO", isByteSwapped ? Py_True: Py_False, sequence_tuple, tuple);
}

static struct PyMethodDef _twoBitIO_methods[] = {
    {"TwoBitIterator",
     (PyCFunction)TwoBitIterator,
     METH_VARARGS | METH_KEYWORDS,
     TwoBitIterator__doc__
    },
    {NULL,          NULL, 0, NULL} /* sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_twoBitIO",
    "Parser for DNA sequence data in 2bit format",
    -1,
    _twoBitIO_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *
PyInit__twoBitIO(void)
{
    PyObject *module;

    if (PyType_Ready(&TwoBitSequenceType) < 0)
        return NULL;

    module = PyModule_Create(&moduledef);
    if (module == NULL) return NULL;

    Py_INCREF(&TwoBitSequenceType);
    PyModule_AddObject(module, "TwoBitSequence", (PyObject*) &TwoBitSequenceType);

    return module;
}
