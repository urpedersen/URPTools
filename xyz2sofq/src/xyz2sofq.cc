//============================================================================
// Name        : xyz2sofq.cpp
// Author      : Ulf R. Pedersen
// Copyright   : GPL, see COPYING for more.
//============================================================================

#include <cstdlib>
#include <cfloat>
/// #include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>

// #include "helpers_io.h"

using namespace std;


//////////// HELPER FUNCTION

// Color conversion from HSV to RGB
// r,g,b values are from 0 to 1
// h = [0,360], s = [0,1], v = [0,1]
//		if s == 0, then h = -1 (undefined)
void hsv2rgb( float *r, float *g, float *b, float h, float s, float v )
{
	int i;
	float f, p, q, t;

	if( s == 0 ) {
		// achromatic (grey)
		*r = *g = *b = v;
		return;
	}

	h /= 60;			// sector 0 to 5
	i = floor( h );
	i %= 6;
	f = h - i;			// factorial part of h
	p = v * ( 1 - s );
	q = v * ( 1 - s * f );
	t = v * ( 1 - s * ( 1 - f ) );

	switch( i ) {
		case 0:
			*r = v;
			*g = t;
			*b = p;
			break;
		case 1:
			*r = q;
			*g = v;
			*b = p;
			break;
		case 2:
			*r = p;
			*g = v;
			*b = t;
			break;
		case 3:
			*r = p;
			*g = q;
			*b = v;
			break;
		case 4:
			*r = t;
			*g = p;
			*b = v;
			break;
		default:		// case 5:
			*r = v;
			*g = p;
			*b = q;
			break;
	}
}


// Convert a number between 0-1 to a color. if <0 then black and >0 then white.
void number2color ( float *hue, float *sat, float *val, float number) {

	// gamma < 1.0 will increes the 
	float gamma = 0.75;
	number = pow(number,gamma);
	
	*hue  = 300-number*495.0;		// Set color hue
	*hue -= 360.0*floor(*hue/360.0);

	if ( number < 0. ) {			// Set color value
		*val=0.;
	} else if ( number < 1./5. ) {			
		*val = sqrt(number*5.);
		//*val = number*5.;
	} else {
		*val = 1.;
	}

	if ( number > 1. ) {		// Set color satuation
		*sat = 0. ;
	} else if (number > 3./4. ) {
		*sat = 1. - ( number - 3./4. ) / ( 1. - 3./4. );
	} else {
		*sat = 1. ;
	}
}







int main ( int argc , char* argv[] )
{


cout << endl << "    ..:: Compute the intermediate scattering S(q) of a xyz-formated configuration ::.. "	<< endl;
cout <<         "                         by Ulf R. Pedersen (2015), http://urp.dk "				<< endl<< endl;

// Initialise
string xyzFile;				// Input coordinated as XYZ file
double boxX=1,boxY=1,boxZ=1;		// Size of periodic box
unsigned long pixelsX=301,pixelsY=301,pixelsZ=301;	// Pixels to use in the final pixture
static double twoPI = 8.0*atan(1.0);


// Handle program options
bool option_COMPUTE_WITH_Z = false;	// Compute all q vectors, i.e. include the z-coordinates = 3D
bool option_CUBIC_BOX = false;     
if( argc == 6 ) {			// 2D mode
	option_COMPUTE_WITH_Z = false;
	option_CUBIC_BOX = false;
	xyzFile = argv[1];
	boxX=atof(argv[2]);
	boxY=atof(argv[3]);
	boxZ=boxX; // dummy value
	pixelsX=atol(argv[4]);
	pixelsY=atol(argv[5]);
	pixelsZ=0;
} else if ( argc == 4 ) {		// 3D mode with a cubic box, this enable binning of S(q)
	option_COMPUTE_WITH_Z = true;
	option_CUBIC_BOX = true;
	xyzFile = argv[1];
	boxX=atof(argv[2]);
	boxY=boxX;
	boxZ=boxX;
	pixelsX=atol(argv[3]);
	pixelsY=pixelsX;
	pixelsZ=pixelsX;
} else if ( argc == 8 ) {		// 3D mode with a orthorompic box
	option_COMPUTE_WITH_Z = true;
        option_CUBIC_BOX = false;
	xyzFile = argv[1];
	boxX=atof(argv[2]);
	boxY=atof(argv[3]);
	boxZ=atof(argv[4]);
	pixelsX=atol(argv[5]);
	pixelsY=atol(argv[6]);
	pixelsZ=atol(argv[7]);
} else {
	cout    << "Usage 2D-mode:              " << argv[0] << " [XYZ file] [X] [Y] [X_pixels] [Y_pixels]" << endl   // TODO get-rid of pixels as input, and use "larges n" as input
		<< "Usage 3D-mode (cubic box):  " << argv[0] << " [XYZ file] [X] [X_pixels] " << endl
		<< "Usage 3D-mode:              " << argv[0] << " [XYZ file] [X] [Y] [Z] [X_pixels] [Y_pixels] [Z_pixels] " << endl
		<< "          * [XYZ file] is a file with the input trajector in the xyz-format." << endl
		<< "          * X, Y and Z gives the size of the periodic box." << endl
		<< "          * X_pixels, Y_pixels and Z_pixels determine largest q-vector computed." << endl 
		<< "          * In 3D-mode, all possible q-vectors are computed," << endl 
		<< "              - a scattering image is not written," << endl
		<< "              - but a Sq_binned.dat file is written if the box is cubic." << endl 
		<< "Notes:"                                         << endl 
		<< "         The type in xyz file must be a integer" << endl 
		<< endl
		<< "Examples:  " << argv[0] << " traj.xyz 20.0 30.0 201 301" << endl
		<< "           " << argv[0] << " traj.xyz 20.0 201" << endl
		<< "           " << argv[0] << " traj.xyz 20.0 19.0 18.0 201 191 181" << endl
		<< endl;
	return 0;
}


// Print informations user about run mode
if ( option_COMPUTE_WITH_Z )
	cout << "Computing all q vectors (3D-mode)." << endl ;
else
	cout << "Computing q-vectors in XY plane (2D-mode)." << endl ;

if ( option_CUBIC_BOX )
	cout << "Periodic box is assumed to be cubic." << endl ;
else
	cout << "Periodic box is assumed to be orthorompic." << endl ;



// Read configuration from xyz input file and compute Sq      TODO Compute Sq with FFT for efficiency
vector<unsigned long> numParticles;	// Vector with number of particles in frames
vector<double> Sq;			// Scatter function
if ( option_COMPUTE_WITH_Z )
{
	for(unsigned int i=0;i<pixelsX*pixelsY*pixelsZ;i++) Sq.push_back(0.0);
} else {
	for(unsigned int i=0;i<pixelsX*pixelsY;i++) Sq.push_back(0.0);
}
long cx=pixelsX/2;			// Center pixel for the scattering
long cy=pixelsY/2;
long cz=pixelsZ/2;
cout << "Pixel with origin: cx= " << cx << " cy= " << cy << " cz= " << cz << endl;
{
	cout << "Computing S(q) on " << xyzFile << ". Frame count:";
	vector<double> x,y,z;			// Particle coordinates
	ifstream ifile;
	ifile.open(xyzFile.c_str(),std::ifstream::in);
	if (!ifile) {
		cerr << endl << "error: " << xyzFile << " could not be opened for reading. Exit." << endl;
		return 1;
	} else {
		string line;
		char *pEnd;

		getline(ifile,line);
		unsigned int frame = 0;
		while( ifile.good() ) {	// Read frames and compute Sq
			if	(frame%50==0 )	cout << endl << frame << "\t|";
			else if (frame%10==0 )	cout << "|";
			else if (frame%5==0 )	cout << ":";
			else if (frame%50==49)	cout << ".|";
			else 			cout << ".";
			cout.flush();

			//cout << "Number of particles line: " << line << endl;
			unsigned long N = atoi(line.c_str());
			numParticles.push_back(N);
			getline(ifile,line);
			//cout << "Comment line: " << line << endl;
			// Read coordinates
			x.clear();y.clear();z.clear();
			for ( unsigned int i = 0 ; i < N ; i++ ) {
				getline(ifile,line);            // TODO  
                                //cout << "Atom line:" << line << endl;
				long this_type = strtol(line.c_str(),&pEnd,10);
				double this_x  = strtod(pEnd,&pEnd);
				double this_y  = strtod(pEnd,&pEnd);
				double this_z  = strtod(pEnd,&pEnd);
                                //cout << "type=" << this_type << " x= " << this_x << " y= " << this_y << " z= " << this_z << endl;
				x.push_back(this_x/boxX);
				y.push_back(this_y/boxY);
				z.push_back(this_z/boxZ);
			}
			
			// Compute Sq	// TODO add openMP to computation loop
			if ( option_COMPUTE_WITH_Z ) // Compute all possible q-vectors
			{
				for(long iz=0; iz<pixelsZ; iz++ )
				{
					for(long iy=0; iy<pixelsY; iy++ )
					{
						for (long ix=0; ix<pixelsX; ix++ )
						{
							double real=0.0, imag=0.0;
							for (unsigned long p=0; p<N; p++ )
							{
								real += cos(twoPI*(ix-cx)*x.at(p)+twoPI*(iy-cy)*y.at(p)+twoPI*(iz-cz)*z.at(p));
								imag -= sin(twoPI*(ix-cx)*x.at(p)+twoPI*(iy-cy)*y.at(p)+twoPI*(iz-cz)*z.at(p));
							}
							Sq.at(iz*pixelsX*pixelsY+iy*pixelsX+ix) += (real*real+imag*imag)/N;
						}
					}
				}
			}
			else
			{
				for(long iy=0; iy<pixelsY; iy++ )
				{
					for (long ix=0; ix<pixelsX; ix++ )
					{
						double real=0.0, imag=0.0;
						for (unsigned long p=0; p<N; p++ )
						{
							real += cos(twoPI*(ix-cx)*x.at(p)+twoPI*(iy-cy)*y.at(p));
							imag -= sin(twoPI*(ix-cx)*x.at(p)+twoPI*(iy-cy)*y.at(p));
						}
						Sq.at(iy*pixelsX+ix) += (real*real+imag*imag)/N;
					}
				}
			}

			frame++;
			getline(ifile,line);
		} // END OF while(ifile.good()) LOOP
		cout << endl;
	}
	ifile.close();

	// Compute the average Sq by deviding by number of frames
	if ( option_COMPUTE_WITH_Z )
	{
		for(unsigned int iz=0; iz<pixelsZ; iz++ )
		{
			for(unsigned int iy=0; iy<pixelsY; iy++ )
			{
				for (unsigned int ix=0; ix<pixelsX; ix++ )
				{
					Sq.at(iz*pixelsX*pixelsY+iy*pixelsX+ix) /= numParticles.size();
				}
			}
		}
	}
	else
	{
		for(unsigned int iy=0; iy<pixelsY; iy++ )
		{
			for (unsigned int ix=0; ix<pixelsX; ix++ )
			{
				Sq.at(iy*pixelsX+ix) /= numParticles.size();
			}
		}
	}
	cout << "Done computing " << numParticles.size() << " frames." << endl;
}

// Do special case of S(k=0) = Variance(N)/N
{
	double Nmean = 0.0;
	for (unsigned long f = 0 ; f < numParticles.size() ; f++ ) 
		Nmean+= numParticles.at(f);
	Nmean/=numParticles.size();
	double N2mean = 0.0;
	for (unsigned long f = 0 ; f < numParticles.size() ; f++ ) 
		N2mean+= numParticles.at(f)*numParticles.at(f);
	N2mean/=numParticles.size();
	double S0 = N2mean/Nmean - Nmean;
	if( option_COMPUTE_WITH_Z )
		Sq.at(cx*pixelsY*pixelsX+cy*pixelsX+cx) = S0;
	else
		Sq.at(cy*pixelsX+cx) = S0;
	cout << "Scattering at origin: S(0)= " << S0 << endl ;
	cout << "   <N>= " << Nmean << " <N**2> = " << N2mean << " std(N)=sqrt(<N**2>-<N>**2)= " << sqrt(N2mean-Nmean*Nmean) << endl;
}




// Write all Sq to ascii data file
ofstream ofile;
ofile.open ("Sq.dat");
if( option_COMPUTE_WITH_Z )
{
	ofile << " # |q| Sq  n_x n_y n_z " << endl;
	for (long iz=0;iz<pixelsZ;iz++)
	{
		for (long iy=0;iy<pixelsY;iy++)
		{
			for (long ix=0;ix<pixelsX;ix++)
			{
				float qx=twoPI*(ix-cx)/boxX;
				float qy=twoPI*(iy-cy)/boxY;
				float qz=twoPI*(iz-cz)/boxZ;
				ofile << sqrt(qx*qx+qy*qy+qz*qz)			<< " " ;
				ofile << Sq.at(iz*pixelsY*pixelsX+iy*pixelsX+ix)	<< " " ;
				ofile << (ix-cx)					<< " " ;
				ofile << (iy-cy)					<< " " ;
				ofile << (iz-cz)					<< " " ;
				ofile << endl ;
			}
		}
	}
}
else
{
	ofile << " # |q| Sq  n_x n_y pixel_x pixel_y" << endl;
	for (long iy=0;iy<pixelsY;iy++)
	{
		for (long ix=0;ix<pixelsX;ix++)
		{
			float qx=twoPI*(ix-cx)/boxX;
			float qy=twoPI*(iy-cy)/boxY;
			ofile << sqrt(qx*qx+qy*qy)	<< " " ;
			ofile << Sq.at(iy*pixelsX+ix)	<< " " ;
			ofile << (ix-cx)		<< " " ;
			ofile << (iy-cy)		<< " " ;
			ofile << ix			<< " " ;
			ofile << iy			<< " " ;
			ofile << endl ;
		}
	}
}
ofile.close();
cout << "Wrote Sq.dat with the computed scattering." << endl;

// Write binned Sq to ascii data file. Can only bin "symmetric" data with below.
if( option_CUBIC_BOX || ( pixelsX==pixelsY && boxX==boxY  ) )
{
	ofstream ofile;
	ofile.open ("Sq_binned.dat");
	ofile << " # |q| <Sq> error std min max numSamples " << endl;

	for( unsigned int r2=0 ; r2<(pixelsX-cx)*(pixelsX-cx)+(pixelsY-cy)*(pixelsY-cy); r2++ )
	{
		unsigned int counter=0;
		double average=0.;
		double std=0.;
		double min=DBL_MAX;
		double max=0.;
		if( option_COMPUTE_WITH_Z )
		{
			for( unsigned int iz=0; iz<pixelsZ; iz++ )
			{
				for( unsigned int iy=0; iy<pixelsY; iy++ )
				{
					for( unsigned int ix=0; ix<pixelsX; ix++ )
					{
						if( r2==(ix-cx)*(ix-cx)+(iy-cy)*(iy-cy)+(iz-cz)*(iz-cz) )
						{
							counter++;
							double this_Sq = Sq.at(iz*pixelsY*pixelsX+iy*pixelsX+ix);
							average+=this_Sq;
							std+=this_Sq*this_Sq;
							if(this_Sq<min)
								min=this_Sq;
							if(Sq.at(iy*pixelsX+ix)>max)
								max=this_Sq;
						}
					}
				}
			}

		}
		else
		{
			for( unsigned int iy=0; iy<pixelsY; iy++ )
			{
				for( unsigned int ix=0; ix<pixelsX; ix++ )
				{
					if( r2==(ix-cx)*(ix-cx)+(iy-cy)*(iy-cy) )
					{
						counter++;
						double this_Sq = Sq.at(iy*pixelsX+ix);
						average+=this_Sq;
						std+=this_Sq*this_Sq;
						if(this_Sq<min)
							min=this_Sq;
						if(this_Sq>max)
							max=this_Sq;
					}
				}
			}
		}
		if( counter>0 )
		{
			double error=0.;
			average/=(double)counter;
			std/=(double)counter;
			std = sqrt(std - average*average)/sqrt(counter);
			if(counter>1)
				error=std/sqrt((double)counter);
			else 
				error=0.;
			ofile << twoPI*sqrt(r2)/boxX		<< " ";
			ofile << average			<< " ";
			ofile << error				<< " ";
			ofile << std				<< " ";
			ofile << min				<< " ";
			ofile << max				<< " ";
			ofile << counter			<< " ";
			ofile << endl ;
		}
	}

	ofile.close();
	cout << "Wrote binned Sq(|q|) to Sq_binned.dat." << endl; 
} else {
	cout	<< "Note: Can only compute the binned Sq(|q|) for in 3D cubic mode. "<< endl ;
}



// Write image in the ppm format
if( !option_COMPUTE_WITH_Z ){
	float SqMaxVal = 12.;
	stringstream out; 
	out << "P3" << endl; 
	out << pixelsX << " " << pixelsY << endl; 
	out << "255" << endl; 
	float red,green,blue; 
	float hue,sat,val; 
	for (unsigned int iy=0;iy<pixelsY;iy++) { 
		for (unsigned int ix=0;ix<pixelsX;ix++) { 			 
			number2color(&hue,&sat,&val, Sq.at(iy*pixelsX+ix)/SqMaxVal );
			hsv2rgb(&red,&green,&blue,hue,sat,val);
			out << (long)floor(red*255) << " " << (long)floor(green*255) << " " << (long)floor(blue*255) << " ";
		}
		out << endl;
	}
	ofstream ofile;
	ofile.open ("Sq.ppm");
	ofile << out.str();
	ofile.close();
	cout << "Wrote 2D scatter spectrum to Sq.ppm with " << SqMaxVal << " as the maximum value on the color scale." << endl;
}
else
{
	cout << "Note: It is not possible to write Sq.ppm image of 3D data." << endl;
}




// Write color bar
if( !option_COMPUTE_WITH_Z ){
	unsigned int colorBarWidth = floor(pixelsY/10+10);
	float red,green,blue; 
	float hue,sat,val; 
	stringstream out; 
	out << "P3" << endl;
	out << colorBarWidth << " " << pixelsY << endl; 
	out << "255" << endl; 
	for (unsigned int iy=0;iy<pixelsY;iy++) { 
		for (unsigned int ix=0;ix<colorBarWidth;ix++) { 
			number2color(&hue,&sat,&val, (double)(pixelsY-iy)/(double)pixelsY );
			hsv2rgb(&red,&green,&blue,hue,sat,val);
			out << (long)floor(red*255) << " " << (long)floor(green*255) << " " << (long)floor(blue*255) << " ";
		}
		out << endl;
	}
	ofstream ofile;
	ofile.open ("ColorBar.ppm");
	ofile << out.str();
	ofile.close();
}




cout << "Happy ending of " << argv[0] << ". Have a nice day!" << endl;
return 0;

}  // END OF MAIN




