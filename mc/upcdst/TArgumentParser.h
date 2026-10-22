// Configuration reader for the interpreted ROOT 5 production macro.
class TArgumentParser: public TMap {
public:
  TArgumentParser() { SetOwner(); }

  void AddBool(const string &name, Bool_t *value) {
    Parameter *par = new Parameter;
    par->boolean = value;
    Add(new TObjString(name.c_str()), par);
  }
  void AddInt(const string &name, Int_t *value) {
    Parameter *par = new Parameter;
    par->integer = value;
    Add(new TObjString(name.c_str()), par);
  }
  void AddFloat(const string &name, Float_t *value) {
    Parameter *par = new Parameter;
    par->floating = value;
    Add(new TObjString(name.c_str()), par);
  }
  void AddFuncInt3(const string &name, void (*value)(UInt_t, Int_t, Int_t)) {
    Parameter *par = new Parameter;
    par->function = value;
    Add(new TObjString(name.c_str()), par);
  }

  void Parse(const string &config) {
    ifstream input(config.c_str());
    string line;
    while(getline(input, line)) {
      string::size_type first = line.find_first_not_of(" \t\r");
      if(first == string::npos || line[first] == '#') continue;
      istringstream values(line);
      string name;
      values >> name;
      TObject *object = GetValue(name.c_str());
      if(!object) continue; // scheduler-only options are not maker parameters
      Parameter *par = static_cast<Parameter*>(object);
      if(par->boolean) values >> *(par->boolean);
      else if(par->integer) values >> *(par->integer);
      else if(par->floating) values >> *(par->floating);
      else if(par->function) {
        UInt_t id = 0;
        Int_t low = 0, high = 0;
        if(!(values >> id)) continue;
        values >> low >> high;
        (*(par->function))(id, low, high);
      }
    }
  }

private:
  // A single concrete wrapper avoids CINT's casts between interpreted template
  // specializations. Exactly one pointer is set when a parameter is registered.
  struct Parameter: public TObject {
    Bool_t *boolean;
    Int_t *integer;
    Float_t *floating;
    void (*function)(UInt_t, Int_t, Int_t);
    Parameter(): boolean(0), integer(0), floating(0) { function = 0; }
  };
};
