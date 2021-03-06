#
# Autogenerated by Thrift Compiler (0.14.0)
#
# DO NOT EDIT UNLESS YOU ARE SURE THAT YOU KNOW WHAT YOU ARE DOING
#
#  options string: py
#

from thrift.Thrift import TType, TMessageType, TFrozenDict, TException, TApplicationException
from thrift.protocol.TProtocol import TProtocolException
from thrift.TRecursive import fix_spec

import sys

from thrift.transport import TTransport
all_structs = []


class SearchMethod(object):
    CV = 1
    CVCustom = 2
    EqualWeighting = 3
    OrderStatistics = 4

    _VALUES_TO_NAMES = {
        1: "CV",
        2: "CVCustom",
        3: "EqualWeighting",
        4: "OrderStatistics",
    }

    _NAMES_TO_VALUES = {
        "CV": 1,
        "CVCustom": 2,
        "EqualWeighting": 3,
        "OrderStatistics": 4,
    }


class DistanceMeasure(object):
    ZScore = 1
    ZScoreHubbinessCorrected = 2
    Correlation = 3

    _VALUES_TO_NAMES = {
        1: "ZScore",
        2: "ZScoreHubbinessCorrected",
        3: "Correlation",
    }

    _NAMES_TO_VALUES = {
        "ZScore": 1,
        "ZScoreHubbinessCorrected": 2,
        "Correlation": 3,
    }


class QueryStatus(object):
    Complete = 1
    Incomplete = 2
    Error = 3

    _VALUES_TO_NAMES = {
        1: "Complete",
        2: "Incomplete",
        3: "Error",
    }

    _NAMES_TO_VALUES = {
        "Complete": 1,
        "Incomplete": 2,
        "Error": 3,
    }


class QueryParams(object):
    """
    Attributes:
     - searchMethod
     - distanceMeasure
     - minQueryGenesFraction
     - minGenomeFraction
     - rbpParam
     - useNegativeCorrelation
     - checkDatasetSize
     - useGeneSymbols
     - simulateWeights

    """


    def __init__(self, searchMethod=1, distanceMeasure=2, minQueryGenesFraction=0.0000000000000000, minGenomeFraction=0.0000000000000000, rbpParam=0.9900000000000000, useNegativeCorrelation=False, checkDatasetSize=False, useGeneSymbols=False, simulateWeights=False,):
        self.searchMethod = searchMethod
        self.distanceMeasure = distanceMeasure
        self.minQueryGenesFraction = minQueryGenesFraction
        self.minGenomeFraction = minGenomeFraction
        self.rbpParam = rbpParam
        self.useNegativeCorrelation = useNegativeCorrelation
        self.checkDatasetSize = checkDatasetSize
        self.useGeneSymbols = useGeneSymbols
        self.simulateWeights = simulateWeights

    def read(self, iprot):
        if iprot._fast_decode is not None and isinstance(iprot.trans, TTransport.CReadableTransport) and self.thrift_spec is not None:
            iprot._fast_decode(self, iprot, [self.__class__, self.thrift_spec])
            return
        iprot.readStructBegin()
        while True:
            (fname, ftype, fid) = iprot.readFieldBegin()
            if ftype == TType.STOP:
                break
            if fid == 1:
                if ftype == TType.I32:
                    self.searchMethod = iprot.readI32()
                else:
                    iprot.skip(ftype)
            elif fid == 2:
                if ftype == TType.I32:
                    self.distanceMeasure = iprot.readI32()
                else:
                    iprot.skip(ftype)
            elif fid == 3:
                if ftype == TType.DOUBLE:
                    self.minQueryGenesFraction = iprot.readDouble()
                else:
                    iprot.skip(ftype)
            elif fid == 4:
                if ftype == TType.DOUBLE:
                    self.minGenomeFraction = iprot.readDouble()
                else:
                    iprot.skip(ftype)
            elif fid == 5:
                if ftype == TType.DOUBLE:
                    self.rbpParam = iprot.readDouble()
                else:
                    iprot.skip(ftype)
            elif fid == 6:
                if ftype == TType.BOOL:
                    self.useNegativeCorrelation = iprot.readBool()
                else:
                    iprot.skip(ftype)
            elif fid == 7:
                if ftype == TType.BOOL:
                    self.checkDatasetSize = iprot.readBool()
                else:
                    iprot.skip(ftype)
            elif fid == 8:
                if ftype == TType.BOOL:
                    self.useGeneSymbols = iprot.readBool()
                else:
                    iprot.skip(ftype)
            elif fid == 9:
                if ftype == TType.BOOL:
                    self.simulateWeights = iprot.readBool()
                else:
                    iprot.skip(ftype)
            else:
                iprot.skip(ftype)
            iprot.readFieldEnd()
        iprot.readStructEnd()

    def write(self, oprot):
        if oprot._fast_encode is not None and self.thrift_spec is not None:
            oprot.trans.write(oprot._fast_encode(self, [self.__class__, self.thrift_spec]))
            return
        oprot.writeStructBegin('QueryParams')
        if self.searchMethod is not None:
            oprot.writeFieldBegin('searchMethod', TType.I32, 1)
            oprot.writeI32(self.searchMethod)
            oprot.writeFieldEnd()
        if self.distanceMeasure is not None:
            oprot.writeFieldBegin('distanceMeasure', TType.I32, 2)
            oprot.writeI32(self.distanceMeasure)
            oprot.writeFieldEnd()
        if self.minQueryGenesFraction is not None:
            oprot.writeFieldBegin('minQueryGenesFraction', TType.DOUBLE, 3)
            oprot.writeDouble(self.minQueryGenesFraction)
            oprot.writeFieldEnd()
        if self.minGenomeFraction is not None:
            oprot.writeFieldBegin('minGenomeFraction', TType.DOUBLE, 4)
            oprot.writeDouble(self.minGenomeFraction)
            oprot.writeFieldEnd()
        if self.rbpParam is not None:
            oprot.writeFieldBegin('rbpParam', TType.DOUBLE, 5)
            oprot.writeDouble(self.rbpParam)
            oprot.writeFieldEnd()
        if self.useNegativeCorrelation is not None:
            oprot.writeFieldBegin('useNegativeCorrelation', TType.BOOL, 6)
            oprot.writeBool(self.useNegativeCorrelation)
            oprot.writeFieldEnd()
        if self.checkDatasetSize is not None:
            oprot.writeFieldBegin('checkDatasetSize', TType.BOOL, 7)
            oprot.writeBool(self.checkDatasetSize)
            oprot.writeFieldEnd()
        if self.useGeneSymbols is not None:
            oprot.writeFieldBegin('useGeneSymbols', TType.BOOL, 8)
            oprot.writeBool(self.useGeneSymbols)
            oprot.writeFieldEnd()
        if self.simulateWeights is not None:
            oprot.writeFieldBegin('simulateWeights', TType.BOOL, 9)
            oprot.writeBool(self.simulateWeights)
            oprot.writeFieldEnd()
        oprot.writeFieldStop()
        oprot.writeStructEnd()

    def validate(self):
        return

    def __repr__(self):
        L = ['%s=%r' % (key, value)
             for key, value in self.__dict__.items()]
        return '%s(%s)' % (self.__class__.__name__, ', '.join(L))

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not (self == other)


class SeekQuery(object):
    """
    Attributes:
     - species
     - genes
     - datasets
     - parameters
     - guideGenes
     - outputDir

    """


    def __init__(self, species="Unknown", genes=None, datasets=None, parameters=None, guideGenes=None, outputDir="/tmp/seek",):
        self.species = species
        self.genes = genes
        self.datasets = datasets
        self.parameters = parameters
        self.guideGenes = guideGenes
        self.outputDir = outputDir

    def read(self, iprot):
        if iprot._fast_decode is not None and isinstance(iprot.trans, TTransport.CReadableTransport) and self.thrift_spec is not None:
            iprot._fast_decode(self, iprot, [self.__class__, self.thrift_spec])
            return
        iprot.readStructBegin()
        while True:
            (fname, ftype, fid) = iprot.readFieldBegin()
            if ftype == TType.STOP:
                break
            if fid == 1:
                if ftype == TType.STRING:
                    self.species = iprot.readString().decode('utf-8', errors='replace') if sys.version_info[0] == 2 else iprot.readString()
                else:
                    iprot.skip(ftype)
            elif fid == 2:
                if ftype == TType.LIST:
                    self.genes = []
                    (_etype3, _size0) = iprot.readListBegin()
                    for _i4 in range(_size0):
                        _elem5 = iprot.readString().decode('utf-8', errors='replace') if sys.version_info[0] == 2 else iprot.readString()
                        self.genes.append(_elem5)
                    iprot.readListEnd()
                else:
                    iprot.skip(ftype)
            elif fid == 3:
                if ftype == TType.LIST:
                    self.datasets = []
                    (_etype9, _size6) = iprot.readListBegin()
                    for _i10 in range(_size6):
                        _elem11 = iprot.readString().decode('utf-8', errors='replace') if sys.version_info[0] == 2 else iprot.readString()
                        self.datasets.append(_elem11)
                    iprot.readListEnd()
                else:
                    iprot.skip(ftype)
            elif fid == 4:
                if ftype == TType.STRUCT:
                    self.parameters = QueryParams()
                    self.parameters.read(iprot)
                else:
                    iprot.skip(ftype)
            elif fid == 5:
                if ftype == TType.LIST:
                    self.guideGenes = []
                    (_etype15, _size12) = iprot.readListBegin()
                    for _i16 in range(_size12):
                        _elem17 = iprot.readString().decode('utf-8', errors='replace') if sys.version_info[0] == 2 else iprot.readString()
                        self.guideGenes.append(_elem17)
                    iprot.readListEnd()
                else:
                    iprot.skip(ftype)
            elif fid == 6:
                if ftype == TType.STRING:
                    self.outputDir = iprot.readString().decode('utf-8', errors='replace') if sys.version_info[0] == 2 else iprot.readString()
                else:
                    iprot.skip(ftype)
            else:
                iprot.skip(ftype)
            iprot.readFieldEnd()
        iprot.readStructEnd()

    def write(self, oprot):
        if oprot._fast_encode is not None and self.thrift_spec is not None:
            oprot.trans.write(oprot._fast_encode(self, [self.__class__, self.thrift_spec]))
            return
        oprot.writeStructBegin('SeekQuery')
        if self.species is not None:
            oprot.writeFieldBegin('species', TType.STRING, 1)
            oprot.writeString(self.species.encode('utf-8') if sys.version_info[0] == 2 else self.species)
            oprot.writeFieldEnd()
        if self.genes is not None:
            oprot.writeFieldBegin('genes', TType.LIST, 2)
            oprot.writeListBegin(TType.STRING, len(self.genes))
            for iter18 in self.genes:
                oprot.writeString(iter18.encode('utf-8') if sys.version_info[0] == 2 else iter18)
            oprot.writeListEnd()
            oprot.writeFieldEnd()
        if self.datasets is not None:
            oprot.writeFieldBegin('datasets', TType.LIST, 3)
            oprot.writeListBegin(TType.STRING, len(self.datasets))
            for iter19 in self.datasets:
                oprot.writeString(iter19.encode('utf-8') if sys.version_info[0] == 2 else iter19)
            oprot.writeListEnd()
            oprot.writeFieldEnd()
        if self.parameters is not None:
            oprot.writeFieldBegin('parameters', TType.STRUCT, 4)
            self.parameters.write(oprot)
            oprot.writeFieldEnd()
        if self.guideGenes is not None:
            oprot.writeFieldBegin('guideGenes', TType.LIST, 5)
            oprot.writeListBegin(TType.STRING, len(self.guideGenes))
            for iter20 in self.guideGenes:
                oprot.writeString(iter20.encode('utf-8') if sys.version_info[0] == 2 else iter20)
            oprot.writeListEnd()
            oprot.writeFieldEnd()
        if self.outputDir is not None:
            oprot.writeFieldBegin('outputDir', TType.STRING, 6)
            oprot.writeString(self.outputDir.encode('utf-8') if sys.version_info[0] == 2 else self.outputDir)
            oprot.writeFieldEnd()
        oprot.writeFieldStop()
        oprot.writeStructEnd()

    def validate(self):
        if self.species is None:
            raise TProtocolException(message='Required field species is unset!')
        if self.genes is None:
            raise TProtocolException(message='Required field genes is unset!')
        return

    def __repr__(self):
        L = ['%s=%r' % (key, value)
             for key, value in self.__dict__.items()]
        return '%s(%s)' % (self.__class__.__name__, ', '.join(L))

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not (self == other)


class StringDoublePair(object):
    """
    Attributes:
     - name
     - value

    """


    def __init__(self, name=None, value=None,):
        self.name = name
        self.value = value

    def read(self, iprot):
        if iprot._fast_decode is not None and isinstance(iprot.trans, TTransport.CReadableTransport) and self.thrift_spec is not None:
            iprot._fast_decode(self, iprot, [self.__class__, self.thrift_spec])
            return
        iprot.readStructBegin()
        while True:
            (fname, ftype, fid) = iprot.readFieldBegin()
            if ftype == TType.STOP:
                break
            if fid == 1:
                if ftype == TType.STRING:
                    self.name = iprot.readString().decode('utf-8', errors='replace') if sys.version_info[0] == 2 else iprot.readString()
                else:
                    iprot.skip(ftype)
            elif fid == 2:
                if ftype == TType.DOUBLE:
                    self.value = iprot.readDouble()
                else:
                    iprot.skip(ftype)
            else:
                iprot.skip(ftype)
            iprot.readFieldEnd()
        iprot.readStructEnd()

    def write(self, oprot):
        if oprot._fast_encode is not None and self.thrift_spec is not None:
            oprot.trans.write(oprot._fast_encode(self, [self.__class__, self.thrift_spec]))
            return
        oprot.writeStructBegin('StringDoublePair')
        if self.name is not None:
            oprot.writeFieldBegin('name', TType.STRING, 1)
            oprot.writeString(self.name.encode('utf-8') if sys.version_info[0] == 2 else self.name)
            oprot.writeFieldEnd()
        if self.value is not None:
            oprot.writeFieldBegin('value', TType.DOUBLE, 2)
            oprot.writeDouble(self.value)
            oprot.writeFieldEnd()
        oprot.writeFieldStop()
        oprot.writeStructEnd()

    def validate(self):
        if self.name is None:
            raise TProtocolException(message='Required field name is unset!')
        if self.value is None:
            raise TProtocolException(message='Required field value is unset!')
        return

    def __repr__(self):
        L = ['%s=%r' % (key, value)
             for key, value in self.__dict__.items()]
        return '%s(%s)' % (self.__class__.__name__, ', '.join(L))

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not (self == other)


class QueryResult(object):
    """
    Attributes:
     - success
     - geneScores
     - datasetWeights
     - status
     - statusMsg

    """


    def __init__(self, success=None, geneScores=None, datasetWeights=None, status=None, statusMsg=None,):
        self.success = success
        self.geneScores = geneScores
        self.datasetWeights = datasetWeights
        self.status = status
        self.statusMsg = statusMsg

    def read(self, iprot):
        if iprot._fast_decode is not None and isinstance(iprot.trans, TTransport.CReadableTransport) and self.thrift_spec is not None:
            iprot._fast_decode(self, iprot, [self.__class__, self.thrift_spec])
            return
        iprot.readStructBegin()
        while True:
            (fname, ftype, fid) = iprot.readFieldBegin()
            if ftype == TType.STOP:
                break
            if fid == 1:
                if ftype == TType.BOOL:
                    self.success = iprot.readBool()
                else:
                    iprot.skip(ftype)
            elif fid == 2:
                if ftype == TType.LIST:
                    self.geneScores = []
                    (_etype24, _size21) = iprot.readListBegin()
                    for _i25 in range(_size21):
                        _elem26 = StringDoublePair()
                        _elem26.read(iprot)
                        self.geneScores.append(_elem26)
                    iprot.readListEnd()
                else:
                    iprot.skip(ftype)
            elif fid == 3:
                if ftype == TType.LIST:
                    self.datasetWeights = []
                    (_etype30, _size27) = iprot.readListBegin()
                    for _i31 in range(_size27):
                        _elem32 = StringDoublePair()
                        _elem32.read(iprot)
                        self.datasetWeights.append(_elem32)
                    iprot.readListEnd()
                else:
                    iprot.skip(ftype)
            elif fid == 4:
                if ftype == TType.I32:
                    self.status = iprot.readI32()
                else:
                    iprot.skip(ftype)
            elif fid == 5:
                if ftype == TType.STRING:
                    self.statusMsg = iprot.readString().decode('utf-8', errors='replace') if sys.version_info[0] == 2 else iprot.readString()
                else:
                    iprot.skip(ftype)
            else:
                iprot.skip(ftype)
            iprot.readFieldEnd()
        iprot.readStructEnd()

    def write(self, oprot):
        if oprot._fast_encode is not None and self.thrift_spec is not None:
            oprot.trans.write(oprot._fast_encode(self, [self.__class__, self.thrift_spec]))
            return
        oprot.writeStructBegin('QueryResult')
        if self.success is not None:
            oprot.writeFieldBegin('success', TType.BOOL, 1)
            oprot.writeBool(self.success)
            oprot.writeFieldEnd()
        if self.geneScores is not None:
            oprot.writeFieldBegin('geneScores', TType.LIST, 2)
            oprot.writeListBegin(TType.STRUCT, len(self.geneScores))
            for iter33 in self.geneScores:
                iter33.write(oprot)
            oprot.writeListEnd()
            oprot.writeFieldEnd()
        if self.datasetWeights is not None:
            oprot.writeFieldBegin('datasetWeights', TType.LIST, 3)
            oprot.writeListBegin(TType.STRUCT, len(self.datasetWeights))
            for iter34 in self.datasetWeights:
                iter34.write(oprot)
            oprot.writeListEnd()
            oprot.writeFieldEnd()
        if self.status is not None:
            oprot.writeFieldBegin('status', TType.I32, 4)
            oprot.writeI32(self.status)
            oprot.writeFieldEnd()
        if self.statusMsg is not None:
            oprot.writeFieldBegin('statusMsg', TType.STRING, 5)
            oprot.writeString(self.statusMsg.encode('utf-8') if sys.version_info[0] == 2 else self.statusMsg)
            oprot.writeFieldEnd()
        oprot.writeFieldStop()
        oprot.writeStructEnd()

    def validate(self):
        if self.success is None:
            raise TProtocolException(message='Required field success is unset!')
        if self.geneScores is None:
            raise TProtocolException(message='Required field geneScores is unset!')
        return

    def __repr__(self):
        L = ['%s=%r' % (key, value)
             for key, value in self.__dict__.items()]
        return '%s(%s)' % (self.__class__.__name__, ', '.join(L))

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not (self == other)
all_structs.append(QueryParams)
QueryParams.thrift_spec = (
    None,  # 0
    (1, TType.I32, 'searchMethod', None, 1, ),  # 1
    (2, TType.I32, 'distanceMeasure', None, 2, ),  # 2
    (3, TType.DOUBLE, 'minQueryGenesFraction', None, 0.0000000000000000, ),  # 3
    (4, TType.DOUBLE, 'minGenomeFraction', None, 0.0000000000000000, ),  # 4
    (5, TType.DOUBLE, 'rbpParam', None, 0.9900000000000000, ),  # 5
    (6, TType.BOOL, 'useNegativeCorrelation', None, False, ),  # 6
    (7, TType.BOOL, 'checkDatasetSize', None, False, ),  # 7
    (8, TType.BOOL, 'useGeneSymbols', None, False, ),  # 8
    (9, TType.BOOL, 'simulateWeights', None, False, ),  # 9
)
all_structs.append(SeekQuery)
SeekQuery.thrift_spec = (
    None,  # 0
    (1, TType.STRING, 'species', 'UTF8', "Unknown", ),  # 1
    (2, TType.LIST, 'genes', (TType.STRING, 'UTF8', False), None, ),  # 2
    (3, TType.LIST, 'datasets', (TType.STRING, 'UTF8', False), None, ),  # 3
    (4, TType.STRUCT, 'parameters', [QueryParams, None], None, ),  # 4
    (5, TType.LIST, 'guideGenes', (TType.STRING, 'UTF8', False), None, ),  # 5
    (6, TType.STRING, 'outputDir', 'UTF8', "/tmp/seek", ),  # 6
)
all_structs.append(StringDoublePair)
StringDoublePair.thrift_spec = (
    None,  # 0
    (1, TType.STRING, 'name', 'UTF8', None, ),  # 1
    (2, TType.DOUBLE, 'value', None, None, ),  # 2
)
all_structs.append(QueryResult)
QueryResult.thrift_spec = (
    None,  # 0
    (1, TType.BOOL, 'success', None, None, ),  # 1
    (2, TType.LIST, 'geneScores', (TType.STRUCT, [StringDoublePair, None], False), None, ),  # 2
    (3, TType.LIST, 'datasetWeights', (TType.STRUCT, [StringDoublePair, None], False), None, ),  # 3
    (4, TType.I32, 'status', None, None, ),  # 4
    (5, TType.STRING, 'statusMsg', 'UTF8', None, ),  # 5
)
fix_spec(all_structs)
del all_structs
