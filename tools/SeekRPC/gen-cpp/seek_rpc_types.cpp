/**
 * Autogenerated by Thrift Compiler (0.14.0)
 *
 * DO NOT EDIT UNLESS YOU ARE SURE THAT YOU KNOW WHAT YOU ARE DOING
 *  @generated
 */
#include "seek_rpc_types.h"

#include <algorithm>
#include <ostream>

#include <thrift/TToString.h>

namespace SeekRPC {

int _kSearchMethodValues[] = {
  SearchMethod::CV,
  SearchMethod::CVCustom,
  SearchMethod::EqualWeighting,
  SearchMethod::OrderStatistics
};
const char* _kSearchMethodNames[] = {
  "CV",
  "CVCustom",
  "EqualWeighting",
  "OrderStatistics"
};
const std::map<int, const char*> _SearchMethod_VALUES_TO_NAMES(::apache::thrift::TEnumIterator(4, _kSearchMethodValues, _kSearchMethodNames), ::apache::thrift::TEnumIterator(-1, nullptr, nullptr));

std::ostream& operator<<(std::ostream& out, const SearchMethod::type& val) {
  std::map<int, const char*>::const_iterator it = _SearchMethod_VALUES_TO_NAMES.find(val);
  if (it != _SearchMethod_VALUES_TO_NAMES.end()) {
    out << it->second;
  } else {
    out << static_cast<int>(val);
  }
  return out;
}

std::string to_string(const SearchMethod::type& val) {
  std::map<int, const char*>::const_iterator it = _SearchMethod_VALUES_TO_NAMES.find(val);
  if (it != _SearchMethod_VALUES_TO_NAMES.end()) {
    return std::string(it->second);
  } else {
    return std::to_string(static_cast<int>(val));
  }
}

int _kDistanceMeasureValues[] = {
  DistanceMeasure::ZScore,
  DistanceMeasure::ZScoreHubbinessCorrected,
  DistanceMeasure::Correlation
};
const char* _kDistanceMeasureNames[] = {
  "ZScore",
  "ZScoreHubbinessCorrected",
  "Correlation"
};
const std::map<int, const char*> _DistanceMeasure_VALUES_TO_NAMES(::apache::thrift::TEnumIterator(3, _kDistanceMeasureValues, _kDistanceMeasureNames), ::apache::thrift::TEnumIterator(-1, nullptr, nullptr));

std::ostream& operator<<(std::ostream& out, const DistanceMeasure::type& val) {
  std::map<int, const char*>::const_iterator it = _DistanceMeasure_VALUES_TO_NAMES.find(val);
  if (it != _DistanceMeasure_VALUES_TO_NAMES.end()) {
    out << it->second;
  } else {
    out << static_cast<int>(val);
  }
  return out;
}

std::string to_string(const DistanceMeasure::type& val) {
  std::map<int, const char*>::const_iterator it = _DistanceMeasure_VALUES_TO_NAMES.find(val);
  if (it != _DistanceMeasure_VALUES_TO_NAMES.end()) {
    return std::string(it->second);
  } else {
    return std::to_string(static_cast<int>(val));
  }
}

int _kQueryStatusValues[] = {
  QueryStatus::Complete,
  QueryStatus::Incomplete,
  QueryStatus::Error
};
const char* _kQueryStatusNames[] = {
  "Complete",
  "Incomplete",
  "Error"
};
const std::map<int, const char*> _QueryStatus_VALUES_TO_NAMES(::apache::thrift::TEnumIterator(3, _kQueryStatusValues, _kQueryStatusNames), ::apache::thrift::TEnumIterator(-1, nullptr, nullptr));

std::ostream& operator<<(std::ostream& out, const QueryStatus::type& val) {
  std::map<int, const char*>::const_iterator it = _QueryStatus_VALUES_TO_NAMES.find(val);
  if (it != _QueryStatus_VALUES_TO_NAMES.end()) {
    out << it->second;
  } else {
    out << static_cast<int>(val);
  }
  return out;
}

std::string to_string(const QueryStatus::type& val) {
  std::map<int, const char*>::const_iterator it = _QueryStatus_VALUES_TO_NAMES.find(val);
  if (it != _QueryStatus_VALUES_TO_NAMES.end()) {
    return std::string(it->second);
  } else {
    return std::to_string(static_cast<int>(val));
  }
}


QueryParams::~QueryParams() noexcept {
}


void QueryParams::__set_search_method(const SearchMethod::type val) {
  this->search_method = val;
__isset.search_method = true;
}

void QueryParams::__set_distance_measure(const DistanceMeasure::type val) {
  this->distance_measure = val;
__isset.distance_measure = true;
}

void QueryParams::__set_min_query_genes_fraction(const double val) {
  this->min_query_genes_fraction = val;
__isset.min_query_genes_fraction = true;
}

void QueryParams::__set_min_genome_fraction(const double val) {
  this->min_genome_fraction = val;
__isset.min_genome_fraction = true;
}

void QueryParams::__set_rbp_param(const double val) {
  this->rbp_param = val;
__isset.rbp_param = true;
}

void QueryParams::__set_useNegativeCorrelation(const bool val) {
  this->useNegativeCorrelation = val;
__isset.useNegativeCorrelation = true;
}

void QueryParams::__set_check_dataset_size(const bool val) {
  this->check_dataset_size = val;
__isset.check_dataset_size = true;
}

void QueryParams::__set_use_gene_symbols(const bool val) {
  this->use_gene_symbols = val;
__isset.use_gene_symbols = true;
}
std::ostream& operator<<(std::ostream& out, const QueryParams& obj)
{
  obj.printTo(out);
  return out;
}


uint32_t QueryParams::read(::apache::thrift::protocol::TProtocol* iprot) {

  ::apache::thrift::protocol::TInputRecursionTracker tracker(*iprot);
  uint32_t xfer = 0;
  std::string fname;
  ::apache::thrift::protocol::TType ftype;
  int16_t fid;

  xfer += iprot->readStructBegin(fname);

  using ::apache::thrift::protocol::TProtocolException;


  while (true)
  {
    xfer += iprot->readFieldBegin(fname, ftype, fid);
    if (ftype == ::apache::thrift::protocol::T_STOP) {
      break;
    }
    switch (fid)
    {
      case 1:
        if (ftype == ::apache::thrift::protocol::T_I32) {
          int32_t ecast0;
          xfer += iprot->readI32(ecast0);
          this->search_method = (SearchMethod::type)ecast0;
          this->__isset.search_method = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 2:
        if (ftype == ::apache::thrift::protocol::T_I32) {
          int32_t ecast1;
          xfer += iprot->readI32(ecast1);
          this->distance_measure = (DistanceMeasure::type)ecast1;
          this->__isset.distance_measure = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 3:
        if (ftype == ::apache::thrift::protocol::T_DOUBLE) {
          xfer += iprot->readDouble(this->min_query_genes_fraction);
          this->__isset.min_query_genes_fraction = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 4:
        if (ftype == ::apache::thrift::protocol::T_DOUBLE) {
          xfer += iprot->readDouble(this->min_genome_fraction);
          this->__isset.min_genome_fraction = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 5:
        if (ftype == ::apache::thrift::protocol::T_DOUBLE) {
          xfer += iprot->readDouble(this->rbp_param);
          this->__isset.rbp_param = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 6:
        if (ftype == ::apache::thrift::protocol::T_BOOL) {
          xfer += iprot->readBool(this->useNegativeCorrelation);
          this->__isset.useNegativeCorrelation = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 7:
        if (ftype == ::apache::thrift::protocol::T_BOOL) {
          xfer += iprot->readBool(this->check_dataset_size);
          this->__isset.check_dataset_size = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 8:
        if (ftype == ::apache::thrift::protocol::T_BOOL) {
          xfer += iprot->readBool(this->use_gene_symbols);
          this->__isset.use_gene_symbols = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      default:
        xfer += iprot->skip(ftype);
        break;
    }
    xfer += iprot->readFieldEnd();
  }

  xfer += iprot->readStructEnd();

  return xfer;
}

uint32_t QueryParams::write(::apache::thrift::protocol::TProtocol* oprot) const {
  uint32_t xfer = 0;
  ::apache::thrift::protocol::TOutputRecursionTracker tracker(*oprot);
  xfer += oprot->writeStructBegin("QueryParams");

  if (this->__isset.search_method) {
    xfer += oprot->writeFieldBegin("search_method", ::apache::thrift::protocol::T_I32, 1);
    xfer += oprot->writeI32((int32_t)this->search_method);
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.distance_measure) {
    xfer += oprot->writeFieldBegin("distance_measure", ::apache::thrift::protocol::T_I32, 2);
    xfer += oprot->writeI32((int32_t)this->distance_measure);
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.min_query_genes_fraction) {
    xfer += oprot->writeFieldBegin("min_query_genes_fraction", ::apache::thrift::protocol::T_DOUBLE, 3);
    xfer += oprot->writeDouble(this->min_query_genes_fraction);
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.min_genome_fraction) {
    xfer += oprot->writeFieldBegin("min_genome_fraction", ::apache::thrift::protocol::T_DOUBLE, 4);
    xfer += oprot->writeDouble(this->min_genome_fraction);
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.rbp_param) {
    xfer += oprot->writeFieldBegin("rbp_param", ::apache::thrift::protocol::T_DOUBLE, 5);
    xfer += oprot->writeDouble(this->rbp_param);
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.useNegativeCorrelation) {
    xfer += oprot->writeFieldBegin("useNegativeCorrelation", ::apache::thrift::protocol::T_BOOL, 6);
    xfer += oprot->writeBool(this->useNegativeCorrelation);
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.check_dataset_size) {
    xfer += oprot->writeFieldBegin("check_dataset_size", ::apache::thrift::protocol::T_BOOL, 7);
    xfer += oprot->writeBool(this->check_dataset_size);
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.use_gene_symbols) {
    xfer += oprot->writeFieldBegin("use_gene_symbols", ::apache::thrift::protocol::T_BOOL, 8);
    xfer += oprot->writeBool(this->use_gene_symbols);
    xfer += oprot->writeFieldEnd();
  }
  xfer += oprot->writeFieldStop();
  xfer += oprot->writeStructEnd();
  return xfer;
}

void swap(QueryParams &a, QueryParams &b) {
  using ::std::swap;
  swap(a.search_method, b.search_method);
  swap(a.distance_measure, b.distance_measure);
  swap(a.min_query_genes_fraction, b.min_query_genes_fraction);
  swap(a.min_genome_fraction, b.min_genome_fraction);
  swap(a.rbp_param, b.rbp_param);
  swap(a.useNegativeCorrelation, b.useNegativeCorrelation);
  swap(a.check_dataset_size, b.check_dataset_size);
  swap(a.use_gene_symbols, b.use_gene_symbols);
  swap(a.__isset, b.__isset);
}

QueryParams::QueryParams(const QueryParams& other2) {
  search_method = other2.search_method;
  distance_measure = other2.distance_measure;
  min_query_genes_fraction = other2.min_query_genes_fraction;
  min_genome_fraction = other2.min_genome_fraction;
  rbp_param = other2.rbp_param;
  useNegativeCorrelation = other2.useNegativeCorrelation;
  check_dataset_size = other2.check_dataset_size;
  use_gene_symbols = other2.use_gene_symbols;
  __isset = other2.__isset;
}
QueryParams& QueryParams::operator=(const QueryParams& other3) {
  search_method = other3.search_method;
  distance_measure = other3.distance_measure;
  min_query_genes_fraction = other3.min_query_genes_fraction;
  min_genome_fraction = other3.min_genome_fraction;
  rbp_param = other3.rbp_param;
  useNegativeCorrelation = other3.useNegativeCorrelation;
  check_dataset_size = other3.check_dataset_size;
  use_gene_symbols = other3.use_gene_symbols;
  __isset = other3.__isset;
  return *this;
}
void QueryParams::printTo(std::ostream& out) const {
  using ::apache::thrift::to_string;
  out << "QueryParams(";
  out << "search_method="; (__isset.search_method ? (out << to_string(search_method)) : (out << "<null>"));
  out << ", " << "distance_measure="; (__isset.distance_measure ? (out << to_string(distance_measure)) : (out << "<null>"));
  out << ", " << "min_query_genes_fraction="; (__isset.min_query_genes_fraction ? (out << to_string(min_query_genes_fraction)) : (out << "<null>"));
  out << ", " << "min_genome_fraction="; (__isset.min_genome_fraction ? (out << to_string(min_genome_fraction)) : (out << "<null>"));
  out << ", " << "rbp_param="; (__isset.rbp_param ? (out << to_string(rbp_param)) : (out << "<null>"));
  out << ", " << "useNegativeCorrelation="; (__isset.useNegativeCorrelation ? (out << to_string(useNegativeCorrelation)) : (out << "<null>"));
  out << ", " << "check_dataset_size="; (__isset.check_dataset_size ? (out << to_string(check_dataset_size)) : (out << "<null>"));
  out << ", " << "use_gene_symbols="; (__isset.use_gene_symbols ? (out << to_string(use_gene_symbols)) : (out << "<null>"));
  out << ")";
}


SeekQuery::~SeekQuery() noexcept {
}


void SeekQuery::__set_species(const std::string& val) {
  this->species = val;
}

void SeekQuery::__set_genes(const std::vector<std::string> & val) {
  this->genes = val;
}

void SeekQuery::__set_datasets(const std::vector<std::string> & val) {
  this->datasets = val;
__isset.datasets = true;
}

void SeekQuery::__set_parameters(const QueryParams& val) {
  this->parameters = val;
__isset.parameters = true;
}

void SeekQuery::__set_guideGenes(const std::vector<std::string> & val) {
  this->guideGenes = val;
__isset.guideGenes = true;
}

void SeekQuery::__set_outputDir(const std::string& val) {
  this->outputDir = val;
__isset.outputDir = true;
}
std::ostream& operator<<(std::ostream& out, const SeekQuery& obj)
{
  obj.printTo(out);
  return out;
}


uint32_t SeekQuery::read(::apache::thrift::protocol::TProtocol* iprot) {

  ::apache::thrift::protocol::TInputRecursionTracker tracker(*iprot);
  uint32_t xfer = 0;
  std::string fname;
  ::apache::thrift::protocol::TType ftype;
  int16_t fid;

  xfer += iprot->readStructBegin(fname);

  using ::apache::thrift::protocol::TProtocolException;

  bool isset_species = false;
  bool isset_genes = false;

  while (true)
  {
    xfer += iprot->readFieldBegin(fname, ftype, fid);
    if (ftype == ::apache::thrift::protocol::T_STOP) {
      break;
    }
    switch (fid)
    {
      case 1:
        if (ftype == ::apache::thrift::protocol::T_STRING) {
          xfer += iprot->readString(this->species);
          isset_species = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 2:
        if (ftype == ::apache::thrift::protocol::T_LIST) {
          {
            this->genes.clear();
            uint32_t _size4;
            ::apache::thrift::protocol::TType _etype7;
            xfer += iprot->readListBegin(_etype7, _size4);
            this->genes.resize(_size4);
            uint32_t _i8;
            for (_i8 = 0; _i8 < _size4; ++_i8)
            {
              xfer += iprot->readString(this->genes[_i8]);
            }
            xfer += iprot->readListEnd();
          }
          isset_genes = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 3:
        if (ftype == ::apache::thrift::protocol::T_LIST) {
          {
            this->datasets.clear();
            uint32_t _size9;
            ::apache::thrift::protocol::TType _etype12;
            xfer += iprot->readListBegin(_etype12, _size9);
            this->datasets.resize(_size9);
            uint32_t _i13;
            for (_i13 = 0; _i13 < _size9; ++_i13)
            {
              xfer += iprot->readString(this->datasets[_i13]);
            }
            xfer += iprot->readListEnd();
          }
          this->__isset.datasets = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 4:
        if (ftype == ::apache::thrift::protocol::T_STRUCT) {
          xfer += this->parameters.read(iprot);
          this->__isset.parameters = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 5:
        if (ftype == ::apache::thrift::protocol::T_LIST) {
          {
            this->guideGenes.clear();
            uint32_t _size14;
            ::apache::thrift::protocol::TType _etype17;
            xfer += iprot->readListBegin(_etype17, _size14);
            this->guideGenes.resize(_size14);
            uint32_t _i18;
            for (_i18 = 0; _i18 < _size14; ++_i18)
            {
              xfer += iprot->readString(this->guideGenes[_i18]);
            }
            xfer += iprot->readListEnd();
          }
          this->__isset.guideGenes = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 6:
        if (ftype == ::apache::thrift::protocol::T_STRING) {
          xfer += iprot->readString(this->outputDir);
          this->__isset.outputDir = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      default:
        xfer += iprot->skip(ftype);
        break;
    }
    xfer += iprot->readFieldEnd();
  }

  xfer += iprot->readStructEnd();

  if (!isset_species)
    throw TProtocolException(TProtocolException::INVALID_DATA);
  if (!isset_genes)
    throw TProtocolException(TProtocolException::INVALID_DATA);
  return xfer;
}

uint32_t SeekQuery::write(::apache::thrift::protocol::TProtocol* oprot) const {
  uint32_t xfer = 0;
  ::apache::thrift::protocol::TOutputRecursionTracker tracker(*oprot);
  xfer += oprot->writeStructBegin("SeekQuery");

  xfer += oprot->writeFieldBegin("species", ::apache::thrift::protocol::T_STRING, 1);
  xfer += oprot->writeString(this->species);
  xfer += oprot->writeFieldEnd();

  xfer += oprot->writeFieldBegin("genes", ::apache::thrift::protocol::T_LIST, 2);
  {
    xfer += oprot->writeListBegin(::apache::thrift::protocol::T_STRING, static_cast<uint32_t>(this->genes.size()));
    std::vector<std::string> ::const_iterator _iter19;
    for (_iter19 = this->genes.begin(); _iter19 != this->genes.end(); ++_iter19)
    {
      xfer += oprot->writeString((*_iter19));
    }
    xfer += oprot->writeListEnd();
  }
  xfer += oprot->writeFieldEnd();

  if (this->__isset.datasets) {
    xfer += oprot->writeFieldBegin("datasets", ::apache::thrift::protocol::T_LIST, 3);
    {
      xfer += oprot->writeListBegin(::apache::thrift::protocol::T_STRING, static_cast<uint32_t>(this->datasets.size()));
      std::vector<std::string> ::const_iterator _iter20;
      for (_iter20 = this->datasets.begin(); _iter20 != this->datasets.end(); ++_iter20)
      {
        xfer += oprot->writeString((*_iter20));
      }
      xfer += oprot->writeListEnd();
    }
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.parameters) {
    xfer += oprot->writeFieldBegin("parameters", ::apache::thrift::protocol::T_STRUCT, 4);
    xfer += this->parameters.write(oprot);
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.guideGenes) {
    xfer += oprot->writeFieldBegin("guideGenes", ::apache::thrift::protocol::T_LIST, 5);
    {
      xfer += oprot->writeListBegin(::apache::thrift::protocol::T_STRING, static_cast<uint32_t>(this->guideGenes.size()));
      std::vector<std::string> ::const_iterator _iter21;
      for (_iter21 = this->guideGenes.begin(); _iter21 != this->guideGenes.end(); ++_iter21)
      {
        xfer += oprot->writeString((*_iter21));
      }
      xfer += oprot->writeListEnd();
    }
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.outputDir) {
    xfer += oprot->writeFieldBegin("outputDir", ::apache::thrift::protocol::T_STRING, 6);
    xfer += oprot->writeString(this->outputDir);
    xfer += oprot->writeFieldEnd();
  }
  xfer += oprot->writeFieldStop();
  xfer += oprot->writeStructEnd();
  return xfer;
}

void swap(SeekQuery &a, SeekQuery &b) {
  using ::std::swap;
  swap(a.species, b.species);
  swap(a.genes, b.genes);
  swap(a.datasets, b.datasets);
  swap(a.parameters, b.parameters);
  swap(a.guideGenes, b.guideGenes);
  swap(a.outputDir, b.outputDir);
  swap(a.__isset, b.__isset);
}

SeekQuery::SeekQuery(const SeekQuery& other22) {
  species = other22.species;
  genes = other22.genes;
  datasets = other22.datasets;
  parameters = other22.parameters;
  guideGenes = other22.guideGenes;
  outputDir = other22.outputDir;
  __isset = other22.__isset;
}
SeekQuery& SeekQuery::operator=(const SeekQuery& other23) {
  species = other23.species;
  genes = other23.genes;
  datasets = other23.datasets;
  parameters = other23.parameters;
  guideGenes = other23.guideGenes;
  outputDir = other23.outputDir;
  __isset = other23.__isset;
  return *this;
}
void SeekQuery::printTo(std::ostream& out) const {
  using ::apache::thrift::to_string;
  out << "SeekQuery(";
  out << "species=" << to_string(species);
  out << ", " << "genes=" << to_string(genes);
  out << ", " << "datasets="; (__isset.datasets ? (out << to_string(datasets)) : (out << "<null>"));
  out << ", " << "parameters="; (__isset.parameters ? (out << to_string(parameters)) : (out << "<null>"));
  out << ", " << "guideGenes="; (__isset.guideGenes ? (out << to_string(guideGenes)) : (out << "<null>"));
  out << ", " << "outputDir="; (__isset.outputDir ? (out << to_string(outputDir)) : (out << "<null>"));
  out << ")";
}


StringDoublePair::~StringDoublePair() noexcept {
}


void StringDoublePair::__set_name(const std::string& val) {
  this->name = val;
}

void StringDoublePair::__set_value(const double val) {
  this->value = val;
}
std::ostream& operator<<(std::ostream& out, const StringDoublePair& obj)
{
  obj.printTo(out);
  return out;
}


uint32_t StringDoublePair::read(::apache::thrift::protocol::TProtocol* iprot) {

  ::apache::thrift::protocol::TInputRecursionTracker tracker(*iprot);
  uint32_t xfer = 0;
  std::string fname;
  ::apache::thrift::protocol::TType ftype;
  int16_t fid;

  xfer += iprot->readStructBegin(fname);

  using ::apache::thrift::protocol::TProtocolException;

  bool isset_name = false;
  bool isset_value = false;

  while (true)
  {
    xfer += iprot->readFieldBegin(fname, ftype, fid);
    if (ftype == ::apache::thrift::protocol::T_STOP) {
      break;
    }
    switch (fid)
    {
      case 1:
        if (ftype == ::apache::thrift::protocol::T_STRING) {
          xfer += iprot->readString(this->name);
          isset_name = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 2:
        if (ftype == ::apache::thrift::protocol::T_DOUBLE) {
          xfer += iprot->readDouble(this->value);
          isset_value = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      default:
        xfer += iprot->skip(ftype);
        break;
    }
    xfer += iprot->readFieldEnd();
  }

  xfer += iprot->readStructEnd();

  if (!isset_name)
    throw TProtocolException(TProtocolException::INVALID_DATA);
  if (!isset_value)
    throw TProtocolException(TProtocolException::INVALID_DATA);
  return xfer;
}

uint32_t StringDoublePair::write(::apache::thrift::protocol::TProtocol* oprot) const {
  uint32_t xfer = 0;
  ::apache::thrift::protocol::TOutputRecursionTracker tracker(*oprot);
  xfer += oprot->writeStructBegin("StringDoublePair");

  xfer += oprot->writeFieldBegin("name", ::apache::thrift::protocol::T_STRING, 1);
  xfer += oprot->writeString(this->name);
  xfer += oprot->writeFieldEnd();

  xfer += oprot->writeFieldBegin("value", ::apache::thrift::protocol::T_DOUBLE, 2);
  xfer += oprot->writeDouble(this->value);
  xfer += oprot->writeFieldEnd();

  xfer += oprot->writeFieldStop();
  xfer += oprot->writeStructEnd();
  return xfer;
}

void swap(StringDoublePair &a, StringDoublePair &b) {
  using ::std::swap;
  swap(a.name, b.name);
  swap(a.value, b.value);
}

StringDoublePair::StringDoublePair(const StringDoublePair& other24) {
  name = other24.name;
  value = other24.value;
}
StringDoublePair& StringDoublePair::operator=(const StringDoublePair& other25) {
  name = other25.name;
  value = other25.value;
  return *this;
}
void StringDoublePair::printTo(std::ostream& out) const {
  using ::apache::thrift::to_string;
  out << "StringDoublePair(";
  out << "name=" << to_string(name);
  out << ", " << "value=" << to_string(value);
  out << ")";
}


QueryResult::~QueryResult() noexcept {
}


void QueryResult::__set_success(const bool val) {
  this->success = val;
}

void QueryResult::__set_gene_scores(const std::vector<StringDoublePair> & val) {
  this->gene_scores = val;
}

void QueryResult::__set_dataset_weights(const std::vector<StringDoublePair> & val) {
  this->dataset_weights = val;
__isset.dataset_weights = true;
}

void QueryResult::__set_status(const QueryStatus::type val) {
  this->status = val;
__isset.status = true;
}

void QueryResult::__set_statusMsg(const std::string& val) {
  this->statusMsg = val;
__isset.statusMsg = true;
}
std::ostream& operator<<(std::ostream& out, const QueryResult& obj)
{
  obj.printTo(out);
  return out;
}


uint32_t QueryResult::read(::apache::thrift::protocol::TProtocol* iprot) {

  ::apache::thrift::protocol::TInputRecursionTracker tracker(*iprot);
  uint32_t xfer = 0;
  std::string fname;
  ::apache::thrift::protocol::TType ftype;
  int16_t fid;

  xfer += iprot->readStructBegin(fname);

  using ::apache::thrift::protocol::TProtocolException;

  bool isset_success = false;
  bool isset_gene_scores = false;

  while (true)
  {
    xfer += iprot->readFieldBegin(fname, ftype, fid);
    if (ftype == ::apache::thrift::protocol::T_STOP) {
      break;
    }
    switch (fid)
    {
      case 1:
        if (ftype == ::apache::thrift::protocol::T_BOOL) {
          xfer += iprot->readBool(this->success);
          isset_success = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 2:
        if (ftype == ::apache::thrift::protocol::T_LIST) {
          {
            this->gene_scores.clear();
            uint32_t _size26;
            ::apache::thrift::protocol::TType _etype29;
            xfer += iprot->readListBegin(_etype29, _size26);
            this->gene_scores.resize(_size26);
            uint32_t _i30;
            for (_i30 = 0; _i30 < _size26; ++_i30)
            {
              xfer += this->gene_scores[_i30].read(iprot);
            }
            xfer += iprot->readListEnd();
          }
          isset_gene_scores = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 3:
        if (ftype == ::apache::thrift::protocol::T_LIST) {
          {
            this->dataset_weights.clear();
            uint32_t _size31;
            ::apache::thrift::protocol::TType _etype34;
            xfer += iprot->readListBegin(_etype34, _size31);
            this->dataset_weights.resize(_size31);
            uint32_t _i35;
            for (_i35 = 0; _i35 < _size31; ++_i35)
            {
              xfer += this->dataset_weights[_i35].read(iprot);
            }
            xfer += iprot->readListEnd();
          }
          this->__isset.dataset_weights = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 4:
        if (ftype == ::apache::thrift::protocol::T_I32) {
          int32_t ecast36;
          xfer += iprot->readI32(ecast36);
          this->status = (QueryStatus::type)ecast36;
          this->__isset.status = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 5:
        if (ftype == ::apache::thrift::protocol::T_STRING) {
          xfer += iprot->readString(this->statusMsg);
          this->__isset.statusMsg = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      default:
        xfer += iprot->skip(ftype);
        break;
    }
    xfer += iprot->readFieldEnd();
  }

  xfer += iprot->readStructEnd();

  if (!isset_success)
    throw TProtocolException(TProtocolException::INVALID_DATA);
  if (!isset_gene_scores)
    throw TProtocolException(TProtocolException::INVALID_DATA);
  return xfer;
}

uint32_t QueryResult::write(::apache::thrift::protocol::TProtocol* oprot) const {
  uint32_t xfer = 0;
  ::apache::thrift::protocol::TOutputRecursionTracker tracker(*oprot);
  xfer += oprot->writeStructBegin("QueryResult");

  xfer += oprot->writeFieldBegin("success", ::apache::thrift::protocol::T_BOOL, 1);
  xfer += oprot->writeBool(this->success);
  xfer += oprot->writeFieldEnd();

  xfer += oprot->writeFieldBegin("gene_scores", ::apache::thrift::protocol::T_LIST, 2);
  {
    xfer += oprot->writeListBegin(::apache::thrift::protocol::T_STRUCT, static_cast<uint32_t>(this->gene_scores.size()));
    std::vector<StringDoublePair> ::const_iterator _iter37;
    for (_iter37 = this->gene_scores.begin(); _iter37 != this->gene_scores.end(); ++_iter37)
    {
      xfer += (*_iter37).write(oprot);
    }
    xfer += oprot->writeListEnd();
  }
  xfer += oprot->writeFieldEnd();

  if (this->__isset.dataset_weights) {
    xfer += oprot->writeFieldBegin("dataset_weights", ::apache::thrift::protocol::T_LIST, 3);
    {
      xfer += oprot->writeListBegin(::apache::thrift::protocol::T_STRUCT, static_cast<uint32_t>(this->dataset_weights.size()));
      std::vector<StringDoublePair> ::const_iterator _iter38;
      for (_iter38 = this->dataset_weights.begin(); _iter38 != this->dataset_weights.end(); ++_iter38)
      {
        xfer += (*_iter38).write(oprot);
      }
      xfer += oprot->writeListEnd();
    }
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.status) {
    xfer += oprot->writeFieldBegin("status", ::apache::thrift::protocol::T_I32, 4);
    xfer += oprot->writeI32((int32_t)this->status);
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.statusMsg) {
    xfer += oprot->writeFieldBegin("statusMsg", ::apache::thrift::protocol::T_STRING, 5);
    xfer += oprot->writeString(this->statusMsg);
    xfer += oprot->writeFieldEnd();
  }
  xfer += oprot->writeFieldStop();
  xfer += oprot->writeStructEnd();
  return xfer;
}

void swap(QueryResult &a, QueryResult &b) {
  using ::std::swap;
  swap(a.success, b.success);
  swap(a.gene_scores, b.gene_scores);
  swap(a.dataset_weights, b.dataset_weights);
  swap(a.status, b.status);
  swap(a.statusMsg, b.statusMsg);
  swap(a.__isset, b.__isset);
}

QueryResult::QueryResult(const QueryResult& other39) {
  success = other39.success;
  gene_scores = other39.gene_scores;
  dataset_weights = other39.dataset_weights;
  status = other39.status;
  statusMsg = other39.statusMsg;
  __isset = other39.__isset;
}
QueryResult& QueryResult::operator=(const QueryResult& other40) {
  success = other40.success;
  gene_scores = other40.gene_scores;
  dataset_weights = other40.dataset_weights;
  status = other40.status;
  statusMsg = other40.statusMsg;
  __isset = other40.__isset;
  return *this;
}
void QueryResult::printTo(std::ostream& out) const {
  using ::apache::thrift::to_string;
  out << "QueryResult(";
  out << "success=" << to_string(success);
  out << ", " << "gene_scores=" << to_string(gene_scores);
  out << ", " << "dataset_weights="; (__isset.dataset_weights ? (out << to_string(dataset_weights)) : (out << "<null>"));
  out << ", " << "status="; (__isset.status ? (out << to_string(status)) : (out << "<null>"));
  out << ", " << "statusMsg="; (__isset.statusMsg ? (out << to_string(statusMsg)) : (out << "<null>"));
  out << ")";
}

} // namespace
