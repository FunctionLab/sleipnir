#ifndef SEEKERROR_H
#define SEEKERROR_H

#include <stdexcept>

#ifdef NDEBUG         // release build
# define DEBUG(msg)
#else                 // debug build
# define DEBUG_MESSAGE(msg)   std::cout << "### DEBUG: " << msg << std::endl;
#endif


#define FILELINE "(" + std::string(__FILE__) + ":" + std::to_string(__LINE__) + "): "


std::string print_exception_stack(const std::exception& e, int level =  0);


// Generic named_error class
class named_error : public std::runtime_error
{
    private:
    std::string name = "Named_Error";
    std::string msg;

    public:
    named_error(std::string err) : std::runtime_error(err) {
        this->msg = name + ": " + err;
    }

    named_error(std::string errName, std::string err) : std::runtime_error(err) {
        this->name = errName;
        this->msg = errName + ": " + err;
    }

    const std::string & errorType() const {
        return this->name;
    }

    virtual const char* what() const noexcept override {
        return this->msg.c_str();
    }
};

// Argument Error
class argument_error : public named_error
{
    public:
    argument_error(std::string err) : named_error("Argument_Error", err) {};
};

// Initialization Error
class init_error : public named_error
{
    public:
    init_error(std::string err) : named_error("Init_Error", err) {};
};

// Configuration Error
class config_error : public named_error
{
    public:
    config_error(std::string err) : named_error("Config_Error", err) {};
};

// Request Error
class request_error : public named_error
{
    public:
    request_error(std::string err) : named_error("Request_Error", err) {};
};

// Query Error
class query_error : public named_error
{
    public:
    query_error(std::string err) : named_error("Query_Error", err) {};
};

// State Error
class state_error : public named_error
{
    public:
    state_error(std::string err) : named_error("State_Error", err) {};
};


#endif  // SEEKERROR_H