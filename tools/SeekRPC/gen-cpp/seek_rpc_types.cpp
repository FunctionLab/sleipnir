/**
 * Autogenerated by Thrift Compiler (0.13.0)
 *
 * DO NOT EDIT UNLESS YOU ARE SURE THAT YOU KNOW WHAT YOU ARE DOING
 *  @generated
 */
#include "seek_rpc_types.h"

#include <algorithm>
#include <ostream>

#include <thrift/TToString.h>

namespace SeekRPC {


QueryParams::~QueryParams() noexcept {
}


void QueryParams::__set_search_method(const std::string& val) {
  this->search_method = val;
__isset.search_method = true;
}

void QueryParams::__set_distance_measure(const std::string& val) {
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
        if (ftype == ::apache::thrift::protocol::T_STRING) {
          xfer += iprot->readString(this->search_method);
          this->__isset.search_method = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 2:
        if (ftype == ::apache::thrift::protocol::T_STRING) {
          xfer += iprot->readString(this->distance_measure);
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
    xfer += oprot->writeFieldBegin("search_method", ::apache::thrift::protocol::T_STRING, 1);
    xfer += oprot->writeString(this->search_method);
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.distance_measure) {
    xfer += oprot->writeFieldBegin("distance_measure", ::apache::thrift::protocol::T_STRING, 2);
    xfer += oprot->writeString(this->distance_measure);
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

QueryParams::QueryParams(const QueryParams& other0) {
  search_method = other0.search_method;
  distance_measure = other0.distance_measure;
  min_query_genes_fraction = other0.min_query_genes_fraction;
  min_genome_fraction = other0.min_genome_fraction;
  rbp_param = other0.rbp_param;
  useNegativeCorrelation = other0.useNegativeCorrelation;
  check_dataset_size = other0.check_dataset_size;
  use_gene_symbols = other0.use_gene_symbols;
  __isset = other0.__isset;
}
QueryParams& QueryParams::operator=(const QueryParams& other1) {
  search_method = other1.search_method;
  distance_measure = other1.distance_measure;
  min_query_genes_fraction = other1.min_query_genes_fraction;
  min_genome_fraction = other1.min_genome_fraction;
  rbp_param = other1.rbp_param;
  useNegativeCorrelation = other1.useNegativeCorrelation;
  check_dataset_size = other1.check_dataset_size;
  use_gene_symbols = other1.use_gene_symbols;
  __isset = other1.__isset;
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
            uint32_t _size2;
            ::apache::thrift::protocol::TType _etype5;
            xfer += iprot->readListBegin(_etype5, _size2);
            this->genes.resize(_size2);
            uint32_t _i6;
            for (_i6 = 0; _i6 < _size2; ++_i6)
            {
              xfer += iprot->readString(this->genes[_i6]);
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
            uint32_t _size7;
            ::apache::thrift::protocol::TType _etype10;
            xfer += iprot->readListBegin(_etype10, _size7);
            this->datasets.resize(_size7);
            uint32_t _i11;
            for (_i11 = 0; _i11 < _size7; ++_i11)
            {
              xfer += iprot->readString(this->datasets[_i11]);
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
            uint32_t _size12;
            ::apache::thrift::protocol::TType _etype15;
            xfer += iprot->readListBegin(_etype15, _size12);
            this->guideGenes.resize(_size12);
            uint32_t _i16;
            for (_i16 = 0; _i16 < _size12; ++_i16)
            {
              xfer += iprot->readString(this->guideGenes[_i16]);
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
    std::vector<std::string> ::const_iterator _iter17;
    for (_iter17 = this->genes.begin(); _iter17 != this->genes.end(); ++_iter17)
    {
      xfer += oprot->writeString((*_iter17));
    }
    xfer += oprot->writeListEnd();
  }
  xfer += oprot->writeFieldEnd();

  if (this->__isset.datasets) {
    xfer += oprot->writeFieldBegin("datasets", ::apache::thrift::protocol::T_LIST, 3);
    {
      xfer += oprot->writeListBegin(::apache::thrift::protocol::T_STRING, static_cast<uint32_t>(this->datasets.size()));
      std::vector<std::string> ::const_iterator _iter18;
      for (_iter18 = this->datasets.begin(); _iter18 != this->datasets.end(); ++_iter18)
      {
        xfer += oprot->writeString((*_iter18));
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
      std::vector<std::string> ::const_iterator _iter19;
      for (_iter19 = this->guideGenes.begin(); _iter19 != this->guideGenes.end(); ++_iter19)
      {
        xfer += oprot->writeString((*_iter19));
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

SeekQuery::SeekQuery(const SeekQuery& other20) {
  species = other20.species;
  genes = other20.genes;
  datasets = other20.datasets;
  parameters = other20.parameters;
  guideGenes = other20.guideGenes;
  outputDir = other20.outputDir;
  __isset = other20.__isset;
}
SeekQuery& SeekQuery::operator=(const SeekQuery& other21) {
  species = other21.species;
  genes = other21.genes;
  datasets = other21.datasets;
  parameters = other21.parameters;
  guideGenes = other21.guideGenes;
  outputDir = other21.outputDir;
  __isset = other21.__isset;
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

StringDoublePair::StringDoublePair(const StringDoublePair& other22) {
  name = other22.name;
  value = other22.value;
}
StringDoublePair& StringDoublePair::operator=(const StringDoublePair& other23) {
  name = other23.name;
  value = other23.value;
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
            uint32_t _size24;
            ::apache::thrift::protocol::TType _etype27;
            xfer += iprot->readListBegin(_etype27, _size24);
            this->gene_scores.resize(_size24);
            uint32_t _i28;
            for (_i28 = 0; _i28 < _size24; ++_i28)
            {
              xfer += this->gene_scores[_i28].read(iprot);
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
            uint32_t _size29;
            ::apache::thrift::protocol::TType _etype32;
            xfer += iprot->readListBegin(_etype32, _size29);
            this->dataset_weights.resize(_size29);
            uint32_t _i33;
            for (_i33 = 0; _i33 < _size29; ++_i33)
            {
              xfer += this->dataset_weights[_i33].read(iprot);
            }
            xfer += iprot->readListEnd();
          }
          this->__isset.dataset_weights = true;
        } else {
          xfer += iprot->skip(ftype);
        }
        break;
      case 4:
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
    std::vector<StringDoublePair> ::const_iterator _iter34;
    for (_iter34 = this->gene_scores.begin(); _iter34 != this->gene_scores.end(); ++_iter34)
    {
      xfer += (*_iter34).write(oprot);
    }
    xfer += oprot->writeListEnd();
  }
  xfer += oprot->writeFieldEnd();

  if (this->__isset.dataset_weights) {
    xfer += oprot->writeFieldBegin("dataset_weights", ::apache::thrift::protocol::T_LIST, 3);
    {
      xfer += oprot->writeListBegin(::apache::thrift::protocol::T_STRUCT, static_cast<uint32_t>(this->dataset_weights.size()));
      std::vector<StringDoublePair> ::const_iterator _iter35;
      for (_iter35 = this->dataset_weights.begin(); _iter35 != this->dataset_weights.end(); ++_iter35)
      {
        xfer += (*_iter35).write(oprot);
      }
      xfer += oprot->writeListEnd();
    }
    xfer += oprot->writeFieldEnd();
  }
  if (this->__isset.statusMsg) {
    xfer += oprot->writeFieldBegin("statusMsg", ::apache::thrift::protocol::T_STRING, 4);
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
  swap(a.statusMsg, b.statusMsg);
  swap(a.__isset, b.__isset);
}

QueryResult::QueryResult(const QueryResult& other36) {
  success = other36.success;
  gene_scores = other36.gene_scores;
  dataset_weights = other36.dataset_weights;
  statusMsg = other36.statusMsg;
  __isset = other36.__isset;
}
QueryResult& QueryResult::operator=(const QueryResult& other37) {
  success = other37.success;
  gene_scores = other37.gene_scores;
  dataset_weights = other37.dataset_weights;
  statusMsg = other37.statusMsg;
  __isset = other37.__isset;
  return *this;
}
void QueryResult::printTo(std::ostream& out) const {
  using ::apache::thrift::to_string;
  out << "QueryResult(";
  out << "success=" << to_string(success);
  out << ", " << "gene_scores=" << to_string(gene_scores);
  out << ", " << "dataset_weights="; (__isset.dataset_weights ? (out << to_string(dataset_weights)) : (out << "<null>"));
  out << ", " << "statusMsg="; (__isset.statusMsg ? (out << to_string(statusMsg)) : (out << "<null>"));
  out << ")";
}

} // namespace