#include <default_gui_model.h>
#include <string>
#include <time.h>

const int Nmax_Episodes = 400; 
const int Nmax_Seeds = 20; 
const int Length_max_pattern = 40000; // i.e. 4s at 10kHz

class SynapticSequences : public DefaultGUIModel
{

  Q_OBJECT

public:
  SynapticSequences(void);
  virtual ~SynapticSequences(void);

  void execute(void);
  void createGUI(DefaultGUIModel::variable_t*, int);
  void customizeGUI(void);

protected:
  virtual void update(DefaultGUIModel::update_flags_t);

private:
  
  void conductance_update(double *, double *);
  double RdmNumber();

  double Vm;
  double Gl;
  double nGe; /* normalized synaptic conductances */
  double nQe;  /* normalized synaptic weight */
  double period;
  double rate;
  int seed, istim;
  double systime;
  long long count;
  double Te, Ee;
  double Fe;
  double Start_Vector[Nmax_Episodes];
  double Stop_Vector[Nmax_Episodes];
  int Seed_Vector[Nmax_Episodes];
  double Patterns[Nmax_Seeds][Length_max_pattern];
  double dt, delay, duration, random_delay;
  int Nseed;
  int step_counter;
  bool stim_on;
  QString root_dir, filename;

  void initParameters();
  void initRandomization(void);
  void storeRandomization(void);
  void set_filename(void);
  double randZeroToOne();
  void LoadConductanceFromFile();
  void fill_conductance_inputs();
};
