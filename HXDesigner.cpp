/*Shell and tube heat exchanger code*/
#include <iostream>
#include <string>
#include <cmath>
#define PI 3.14159265

using namespace std;

struct shellSide{
	double pressure;
	double designPressure;
	int th;
	int tmin;
	double F;
	double J;
	int ID;
	int OD;
	double ca;
};

struct head{
	double pressure;
	double designPressure;
	int th;
	int tmin;
	double F;
	double J;
	int Rc;
	int Rk;
	double ca;
	float W;
};

struct tubeSide{
	int number;
	int passes;
	double pressure;
	double designPressure;
	int th;
	int tmin;
	double F;
	int ID;
	int OD;
	double ca;
	long E;
	int pitch;
	string pitchType;
	float mu;
	int Pc;
	int l;
};

struct gasket{
	int ID;
	int OD;
	int G;
	float m;
	double Ya;
	double b0;
	double b;
	int bcd;
	float k;
	long H;
	int h_g;
};

struct baffle{
	int cut,length;
};

struct bolt{
	int number;
	int size;
	long W_m1;
	long W_m2;
	long W_m;
	int F;
	double A_m;
	double A;
	int d;
	long h;
};

struct flange{
	int th;
};

struct channel{
	string gtype;
	float kc;
	int th;
};

int main(){
	struct shellSide shell;
	struct tubeSide tube;
	cout<<"Enter the shell and tube material corrosion allowance:"<<endl;
	cin>>shell.ca>>tube.ca;
	cout<<endl;
	cout<<"Enter the number of tubes:"<<endl;
	cin>>tube.number;
	cout<<endl;
	cout<<"Enter the shell side and tube side pressures:"<<endl;
	cin>>shell.pressure>>tube.pressure;
	cout<<endl;
	shell.designPressure=1.1*shell.pressure;
	tube.designPressure=1.1*tube.pressure;
	
	cout<<endl;
	cout<<"Enter the tube pitch:"<<endl;
	cin>>tube.pitch;
	cout<<endl;
	
	cout<<"Enter the external diameter of the tube:"<<endl;
	cin>>tube.OD;
	cout<<endl;
	cout<<"Enter the number of tubes and passes on the tube side:"<<endl;
	cin>>tube.number>>tube.passes;
	cout<<endl;
	cout<<"Enter the type of pitch used in the exchanger:"<<endl;
	cin>>tube.pitchType;
	shellIDprocedure:
		if(tube.pitchType.compare("Square")==0||tube.pitchType.compare("square")==0||tube.pitchType.compare("SQUARE")==0){ //Internal diameter for aquare pitch
			shell.ID=2*tube.pitch*sqrt(tube.passes*tube.number/PI);
		}else if(tube.pitchType.compare("Triangular")==0||tube.pitchType.compare("triangular")==0||tube.pitchType.compare("TRIANGULAR")==0){ //Internal diameter for triangular pitch
			shell.ID=tube.pitch*sqrt(2*tube.passes*tube.number*sqrt(3)/PI);
		}else{
			goto shellIDprocedure; //The pitch type entered is invalid so program reverts to shell ID determination
		}
	
	shell.th=static_cast<int>(shell.designPressure*shell.ID/((2*shell.F*shell.J)-shell.designPressure));
	if(shell.th>shell.tmin)
		shell.th=shell.th;
	else
		shell.th=shell.tmin;
	shell.th=shell.th+shell.ca;
	
	struct head HXhead;
	HXhead.Rc=shell.ID;
	cout<<"Enter the knuckle radius:"<<endl;
	cin>>HXhead.Rk;
	WfactorComp:
		if(HXhead.Rk==0){
			cerr<<"Invalid mathematical calculation\n";
			goto WfactorComp;
		}else
			HXhead.W=0.25*(3+sqrt(HXhead.Rc/HXhead.Rk));
		
	HXhead.th=static_cast<int>(shell.designPressure*HXhead.Rc/(2*shell.F*shell.J))+shell.ca;
	
	struct baffle bfl;
	cout<<"Enter the desired baffle cut for the heat exchanger:"<<endl;
	cin>>bfl.cut;
	bfl.length=static_cast<int>(shell.ID*(1-(0.01*static_cast<float>(bfl.cut))));
	
	struct gasket G1;
	cout<<"Enter the gasket internal diameter:"<<endl;
	cin>>G1.ID;
	cout<<"Enter the gasket external diameter:"<<endl;
	cin>>G1.OD;
	cout<<"Enter the gasket factor:"<<endl;
	cin>>G1.m;
	cout<<"Enter the gasket seating stress:"<<endl;
	cin>>G1.Ya;
	
	G1.G=(int)((G1.ID+G1.OD)/2);
	G1.b0=(G1.OD-G1.ID)/4;
	if(G1.b0>6.3){
		G1.b0=2.5*sqrt(G1.b0);
	}else{
		G1.b0=G1.b0;
	}
	
	struct bolt b1;
	b1.number=static_cast<int>(G1.G/25);
	if(b1.number%4!=0){
		switch(b1.number%4){
			case 1:
				b1.number-=1;
				break;
			case 2:
				b1.number+=2;
				break;
			case 3:
				b1.number+=1;
				break;
			default:
				break;
		}
	}else{
		b1.number=b1.number;
	}
	
	cout<<"Enter the maximum bolt load:"<<endl;
	cin>>b1.F;
	cout<<endl;
	b1.W_m1=PI*G1.b*G1.G*G1.Ya;
	b1.W_m2=(2*PI*G1.b*G1.m*shell.designPressure)+(PI*pow(G1.G,2)*shell.designPressure/4);
	b1.W_m=(b1.W_m1>b1.W_m2)?b1.W_m1:b1.W_m2;
	b1.A_m=b1.W_m/b1.F;
	b1.A=b1.A_m/b1.number;
	b1.d=static_cast<int>((1/0.51)*pow(b1.A,1/2.09));
	
	G1.bcd=G1.G+(2*b1.d)+(2*G1.b0);
	G1.H=PI*pow(G1.G,2)*shell.designPressure/4;
	G1.h_g=(G1.bcd-G1.G)/2;
	G1.k=1/(0.3+(1.5*b1.W_m*G1.h_g/(G1.H*G1.G)));
	
	struct flange fl;
	fl.th=G1.G*sqrt(shell.designPressure/(G1.k*shell.F));
	
	cout<<"Enter the external diameter of tube:"<<endl;
	cin>>tube.OD;
	cout<<"Enter the ultimate stress of tube material:"<<endl;
	cin>>tube.F;
	cout<<"Enter the elasticity modulus of tube material:"<<endl;
	cin>>tube.E;
	cout<<"Enter the Poisson ratio of tube material:"<<endl;
	cin>>tube.mu;
	cout<<"Enter the length of tube:"<<endl;
	cin>>tube.l;
	tube.th=static_cast<int>(tube.designPressure*tube.OD/((2*tube.F)-tube.designPressure));
	tubePresCorrection:
		tube.Pc=static_cast<int>(2.42*tube.E*pow(tube.th/tube.OD,2.5)/(pow(1-pow(tube.mu,2),0.75)*((tube.l/tube.OD)-0.45*sqrt(tube.th/tube.OD))));
		tube.Pc=tube.Pc/4;
		if(shell.designPressure>tube.Pc){
			tube.th+=0.1;
			goto tubePresCorrection;
		}
	
	string exchanger_type;
	cout<<"Enter the type of tube sheet, can be fixed tube sheet or U-tube type:"<<endl;
	cin>>exchanger_type;
	float F;
	if(exchanger_type.compare("fixed")==0||exchanger_type.compare("Fixed")==0||exchanger_type.compare("FIXED")==0){
		F=1;
	}else if(exchanger_type.compare("u-tube")==0||exchanger_type.compare("U-tube")==0||exchanger_type.compare("U-TUBE")==0){
		F=1.25;
	}
	float t_tubesheet=F*G1.G*sqrt(0.25*tube.designPressure/tube.F);
	
	struct channel ch;
	cout<<"Enter the type of gasket being used, can be full facing or ring type:"<<endl;
	cin>>ch.gtype;
	if(ch.gtype.compare("ring")==0||ch.gtype.compare("Ring")==0||ch.gtype.compare("RING")==0){
		ch.kc=0.25;
	}else if(ch.gtype.compare("full-facing")==0||ch.gtype.compare("Full-facing")==0||ch.gtype.compare("FULL-FACING")==0){
		ch.kc=0.3;
	}
	ch.th=static_cast<int>(G1.G*sqrt(ch.kc*tube.designPressure/tube.F));
	
	struct gasket G2;
	struct bolt b2;
	G2.G=G1.G;
	G2.m=G1.m;
	G2.Ya=G1.Ya;
	G2.b0=G1.b0/4;
	b2.number=G2.G/25;
	if(b2.number%4!=0){
		switch(b2.number%4){
			case 1:
				b2.number=b2.number-1;
				break;
			case 2:
				b2.number=b2.number+2;
			    break;
			case 3:
				b2.number=b2.number+1;
				break;
		}
	}else{
		b2.number=b2.number;
	}
	
	b2.W_m1=PI*G2.b*G1.G*G2.Ya;
	b2.W_m2=(2*PI*G2.b*G1.m*shell.designPressure)+(PI*pow(G2.G,2)*shell.designPressure/4);
	b2.W_m=(b2.W_m1>b2.W_m2)?b2.W_m1:b2.W_m2;
	b2.A_m=b2.W_m/b2.F;
	b2.A=b2.A_m/b2.number;
	b2.d=static_cast<int>((1/0.51)*pow(b2.A,1/2.09));
	
	G2.bcd=G2.G+(2*b2.d)+(2*G2.b0);
	G2.H=PI*pow(G2.G,2)*shell.designPressure/4;
	G2.h_g=(G2.bcd-G2.G)/2;
	G2.k=1/(0.3+(1.5*b2.W_m*G2.h_g/(G2.H*G2.G)));
	
	struct flange f2;
	f2.th=static_cast<int>(G2.G*sqrt(shell.designPressure/(G2.k*shell.F)));
	
	int nozzleID;
	cout<<"Enter the internal diameter of nozzle:"<<endl;
	cin>>nozzleID;
	int tnozzle=static_cast<int>(tube.designPressure*nozzleID/((2*tube.F)-tube.designPressure));
	int tionozzle=static_cast<int>(shell.designPressure*nozzleID/((2*shell.F*shell.J)-shell.designPressure)+shell.ca);
	
	cout<<"Unless mentioned as a number, all the computed variables are in mm"<<endl;
	cout<<"The shell ID is "<<shell.ID<<endl;
	cout<<"The shell thickness is :"<<shell.th<<endl;
	cout<<"The thickness of inlet-outlet nozzles is :"<<tionozzle<<endl;
	cout<<"The head is of thickness "<<HXhead.th<<" despite the heuristic measure of "<<shell.th<<endl;
	cout<<"The shell and tube-sheet flange joint has a thickness of "<<fl.th<<endl;
	cout<<"The tube thickness is: "<<tube.th<<endl;
	cout<<"The thickness of tubesheet is: "<<t_tubesheet<<endl;
	cout<<"The thickness of channel is: "<<ch.th<<endl;
	cout<<"The tube sheet and channel flange thickness is:"<<f2.th<<endl;
	int bolts=b1.number>b2.number?b1.number:b2.number;
	cout<<"The number of bolts required is: "<<bolts<<endl;
	return 0;
}