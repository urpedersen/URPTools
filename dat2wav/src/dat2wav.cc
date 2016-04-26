/**
 *  data2wav - convert data stream into an audio file in the uncompressed PCM WAV format
 *  
 *     v0.2    2013-10-08   Ulf R. Pedersen (http://urp.dk)
 *
 *  Usage help: 
 *     data2wav -h
 * 
*/

#include <fstream> 
#include <vector> 
#include <iostream> 
#include <string> 
#include <sstream> 
#include <climits> 
#include <stdlib.h>

using namespace std;

int main(int argc,char* argv[]) { 

   // Default values
   string ofile = "out.wav";
   int SampleRate = 44100;
   int userBitRate = 16;
   short NumChannels = 1;

   // Read command line input
   for (int i = 1; i < argc; i++) {
      if(string(argv[i])=="-h" || string(argv[i])=="--help"){
	      cout << endl 
		      << "   dat2wav v0.2 - convert data stream into an audio file in the uncompressed PCM WAV format (out.wav)" << endl <<  endl
		      << "Usage:"<< endl << endl
		      << "   The program expects a standard input with one float value on each line. The stream should be ended with EOF (ctrl-d). The input data is rescaled so the volume of the wav-output is adjusted to the range of the input. Optional flags are:" << endl << endl
          << " -h or --help" << endl
          << "     Print this help message." << endl << endl
		      << " -o [output filename]" << endl
          << "     Default: out.wav" << endl << endl
		      << " -r [sample rate in Hz]" << endl
          << "     Default: 44100. Recommended values are 96000, 48000, 44100, 32000, 22050, 16000, 11025 or 8000." << endl << endl
		      << " -b [bits per sample]" << endl
          << "     Default: 16. Accepted values are 8, 16 or 32." << endl << endl
		      << " -c [number of channels]" << endl
          << "     Default: 1. Use 2 for stereo. Stereo input is given as L1 R1 L2 R2 L3 R3 ... etc." << endl
		      << endl << endl
		      << "Usage examples: " << endl << endl
                      << "   echo -e \"0.1\\n0.2\\n0.3\" | " << argv[0] << endl << endl
		      << "   cat myData.ascii | " << argv[0] << " -o myData.wav" << endl << endl
                      << "   paste left.dat right.dat | awk '{print $1;print $2}'  | " << argv[0] << " -c 2 -o stereo.wav" << endl << endl
                      << "   awk -F, '{print $2}' data.csv | " << argv[0] << endl << endl
		      << "   awk 'BEGIN{twopi=8.*atan2(1.,1.);for(n=0;n<48000;n++){print sin(twopi*n/48000*110*2^(7/12))}}' | " << argv[0] << " -r 48000 -o E3.wav" << endl << endl
		      << "   awk 'BEGIN{twopi=8.*atan2(1.,1.);for(n=0;n<32000;n++){print sin(twopi*n/32000*110*2^(7/12));print sin(twopi*n/32000*220*2^(7/12))}}' | " << argv[0] << " -r 32000 -c 2 -o E3E4_stereo.wav" << endl << endl << endl
		      << "Copyright (C) 2013 Ulf R. Pedersen (http://urp.dk). "
		      << "This program comes with ABSOLUTELY NO WARRANTY. "
		      << "This is open-source free software (GNU GPL 3.0), and you are welcome to redistribute it "
		      << "under certain conditions (please read LICENSE details)." << endl << endl 
	        ;exit (0);
      }
      else if(string(argv[i])=="-o" && (argc-1)!=i){ // Set name of output file
           i++;
           ofile = string(argv[i]);
      }
      else if(string(argv[i])=="-r" && (argc-1)!=i){ // Set sample rate 
           i++;
           SampleRate = atoi(argv[i]);
      }
      else if(string(argv[i])=="-c" && (argc-1)!=i){ // Set number of channels
           i++;
           NumChannels = atoi(argv[i]);
      }
      else if(string(argv[i])=="-b" && (argc-1)!=i){ // Set number of channels
           i++;
           userBitRate = atoi(argv[i]);
           if(!( userBitRate==8 || userBitRate==16 || userBitRate==32 )){
              cout<<"Error in input bitrate (-b). Use -h for help."<<endl;exit (1);
           }
      }
      else{
         cout<<"Error in input arguments. Use -h for help."<<endl;exit (1);
      }
   }




   // Read data from standard input and compute value_avg and value_del for later rescaling of data
   vector<float> data;
   float value=0.0;
   float value_avg=0;
   string tmpstr;
   while (!cin.eof())
   {
        getline (cin,tmpstr);
        stringstream(tmpstr) >> value;
        //cout << value << endl;
        data.push_back(value);
        value_avg+=value;
   }
   value_avg/=(float)data.size();
   //cout << "avg = " << value_avg << endl;

   float value_min=value_avg; 
   float value_max=value_avg; 
   for ( std::vector<float>::iterator it = data.begin(); it != data.end(); ++it ) {
      if(*it>value_max) value_max=*it;
      if(*it<value_min) value_min=*it;
   }
   //cout << "min = " << value_min << endl;
   //cout << "max = " << value_max << endl;
   float value_del = value_max-value_avg;
   if( value_avg-value_min>value_del ) value_del = value_avg-value_min;
   //cout << "del = " << value_del << endl;


  float max_volume = 0.90;  // [0.;1.] determining the maximum volume 


if(userBitRate==32){
   
   // Fill integer buffer for WAV file with scaled values of the float standard input
   vector<signed int> buffer; 
   int BytesPerSample = sizeof(signed int); // signed int = 4 byte = 32 bit
   int BitsPerSample = BytesPerSample*8;

   signed int ival;
   for ( std::vector<float>::iterator it = data.begin(); it != data.end(); ++it ) {
       ival = (signed int)((*it-value_avg)*max_volume/value_del*INT_MAX);
       //cout << ival << endl;
       buffer.push_back(ival);
   }
   
   int NumSamples = buffer.size();
   // cout << "NumSamples = " << NumSamples << endl;

   // Write the WAV file. Documentation: https://ccrma.stanford.edu/courses/422/projects/WaveFormat/
   std::ofstream stream(ofile.c_str(), std::ios::binary);

   stream.write("RIFF", 4);   // header chunk
   int ChunkSize=(36+BytesPerSample*NumSamples);
   stream.write((const char*)&ChunkSize, 4); 
   stream.write("WAVE", 4);

   stream.write("fmt ", 4);   // formating chunk
   int Subchunk1Size = 16;
   stream.write((const char*)&Subchunk1Size, 4);
   short AudioFormat = 1;     // 1 = PCM 
   stream.write((const char*)&AudioFormat, 2);
   stream.write((const char*)&NumChannels, 2);
   stream.write((const char*)&SampleRate, 4);
   int ByteRate = SampleRate * NumChannels * BytesPerSample; // SampleRate * NumChannels * BytesPerSample
   stream.write((const char*)&ByteRate, 4);
   short BlockAlign = NumChannels * BitsPerSample; // 
   stream.write((const char*)&BlockAlign, 2);
   stream.write((const char*)&BitsPerSample, 2);

   stream.write("data", 4);   // data chunk
   int Subchunk2Size = NumSamples * BytesPerSample; // NumSamples * NumChannels * BytesPerSample
   stream.write((const char*)&Subchunk2Size, 4);
   stream.write((const char*)&buffer.at(0),Subchunk2Size);

   stream.close();


} else if (userBitRate==16) { 

   //cout << "16 bit" << endl;

  // Fill integer buffer for WAV file with scaled values of the float standard input
   vector<signed short> buffer; 
   int BytesPerSample = sizeof(signed short); // signed short = 2 byte = 16 bit
   int BitsPerSample = BytesPerSample*8;

   //cout << BitsPerSample << endl;

   signed short ival;
   for ( std::vector<float>::iterator it = data.begin(); it != data.end(); ++it ) {
       ival = (signed short)((*it-value_avg)*max_volume/value_del*SHRT_MAX);
       //cout << ival << endl;
       buffer.push_back(ival);
   }
   
   int NumSamples = buffer.size();
   // cout << "NumSamples = " << NumSamples << endl;

   // Write the WAV file. Documentation: https://ccrma.stanford.edu/courses/422/projects/WaveFormat/
   std::ofstream stream(ofile.c_str(), std::ios::binary);

   stream.write("RIFF", 4);   // header chunk
   int ChunkSize=(36+BytesPerSample*NumSamples);
   stream.write((const char*)&ChunkSize, 4); 
   stream.write("WAVE", 4);

   stream.write("fmt ", 4);   // formating chunk
   int Subchunk1Size = 16;
   stream.write((const char*)&Subchunk1Size, 4);
   short AudioFormat = 1;     // 1 = PCM 
   stream.write((const char*)&AudioFormat, 2);
   stream.write((const char*)&NumChannels, 2);
   stream.write((const char*)&SampleRate, 4);
   int ByteRate = SampleRate * NumChannels * BytesPerSample; // SampleRate * NumChannels * BytesPerSample
   stream.write((const char*)&ByteRate, 4);
   short BlockAlign = NumChannels * BitsPerSample; // 
   stream.write((const char*)&BlockAlign, 2);
   stream.write((const char*)&BitsPerSample, 2);

   stream.write("data", 4);   // data chunk
   int Subchunk2Size = NumSamples * BytesPerSample; // NumSamples * NumChannels * BytesPerSample
   stream.write((const char*)&Subchunk2Size, 4);
   stream.write((const char*)&buffer.at(0),Subchunk2Size);

   stream.close();


} else if (userBitRate==8) {   // Note. The wave format breaks it's convention for 8bit, and use an unsigned integrer.

   //cout << "8 bit" << endl;

  // Fill integer buffer for WAV file with scaled values of the float standard input
   vector<unsigned char> buffer; 
   int BytesPerSample = sizeof(unsigned char); // signed char = 1 byte = 8 bit
   int BitsPerSample = BytesPerSample*8;

   //cout << BitsPerSample << endl;

   unsigned char ival;
   for ( std::vector<float>::iterator it = data.begin(); it != data.end(); ++it ) {
       ival = (unsigned char)((*it-value_avg+value_del+(1.-max_volume)*value_del)*max_volume/value_del*.5*UCHAR_MAX);
       //cout << ival << endl;
       buffer.push_back(ival);
   }
   
   int NumSamples = buffer.size();
   // cout << "NumSamples = " << NumSamples << endl;

   // Write the WAV file. Documentation: https://ccrma.stanford.edu/courses/422/projects/WaveFormat/
   std::ofstream stream(ofile.c_str(), std::ios::binary);

   stream.write("RIFF", 4);   // header chunk
   int ChunkSize=(36+BytesPerSample*NumSamples);
   stream.write((const char*)&ChunkSize, 4); 
   stream.write("WAVE", 4);

   stream.write("fmt ", 4);   // formating chunk
   int Subchunk1Size = 16;
   stream.write((const char*)&Subchunk1Size, 4);
   short AudioFormat = 1;     // 1 = PCM 
   stream.write((const char*)&AudioFormat, 2);
   stream.write((const char*)&NumChannels, 2);
   stream.write((const char*)&SampleRate, 4);
   int ByteRate = SampleRate * NumChannels * BytesPerSample; // SampleRate * NumChannels * BytesPerSample
   stream.write((const char*)&ByteRate, 4);
   short BlockAlign = NumChannels * BitsPerSample; // 
   stream.write((const char*)&BlockAlign, 2);
   stream.write((const char*)&BitsPerSample, 2);

   stream.write("data", 4);   // data chunk
   int Subchunk2Size = NumSamples * BytesPerSample; // NumSamples * NumChannels * BytesPerSample
   stream.write((const char*)&Subchunk2Size, 4);
   stream.write((const char*)&buffer.at(0),Subchunk2Size);

   stream.close();



} else {
   cout<<"Error. Unknown bitrate. Use -h for help."<<endl;exit (1);
}

   return 0;
}


