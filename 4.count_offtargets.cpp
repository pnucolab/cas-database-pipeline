#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <condition_variable>

#define FILE_CNT 21
#define MAX_THREAD_CNT 64

using namespace std;

class semaphore
{
private:
    mutex mutex_;
    condition_variable condition_;
    unsigned long count_;

public:
    semaphore(int count)
        : count_(0)
    {
        count_ = count;
    }

    void notify()
    {
        unique_lock<mutex> lock(mutex_);
        ++count_;
        condition_.notify_one();
    }

    void wait()
    {
        unique_lock<mutex> lock(mutex_);
        while(!count_)
            condition_.wait(lock);
        --count_;
    }
};

#ifndef stoi
int stoi(string &s) {
    return atoi(s.c_str());
}
#endif

vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

class Info_Transcript {
    public:
        int pk;
        vector< pair<int, int> > cds_pos;
    
        Info_Transcript(int arg_pk, vector<string> cds_strs) {
           int i;
           for (i=0; i<cds_strs.size(); i++) {
               vector<string> splitted_cds_str = split(cds_strs[i], ',');
               cds_pos.push_back(make_pair(stoi(splitted_cds_str[0]), stoi(splitted_cds_str[1])));
           }
           pk = arg_pk;
        }
};

int find_transcript(vector<Info_Transcript*>* transcripts, int pos, char direction, int seqlen=23) {
    int i, j, bp;
    if (direction == '+')
        bp = pos + seqlen - 6;
    else
        bp = pos + 6;
    for (i=0; i<transcripts->size(); i++) {
        for (j=0; j<(*transcripts)[i]->cds_pos.size(); j++) {
            if ((*transcripts)[i]->cds_pos[j].first < bp && bp < (*transcripts)[i]->cds_pos[j].second) {
                return (*transcripts)[i]->pk;
            }
        }
    }
    return -1;
}

static inline string &rstrip(string &s, char delim=0) {
    if ( (delim > 0 && s[s.size()-1] == delim) || (delim == 0 && ( s[s.size()-1] == '\n' ||
                                                                   s[s.size()-1] == '\r' ||
                                                                   s[s.size()-1] == '\t' ||
                                                                   s[s.size()-1] == ' ' ) ) ) {
        s.erase(s.size()-1);
        return rstrip(s, delim);
    }
    return s;
}

vector<string> slice(const vector<string>& v, int start=0, int end=-1) {
    int oldlen = v.size();
    int newlen;

    if (end == -1 or end >= oldlen){
        newlen = oldlen-start;
    } else {
        newlen = end-start;
    }

    vector<string> nv(newlen);

    for (int i=0; i<newlen; i++) {
        nv[i] = v[start+i];
    }
    return nv;
}

class OffCount {
    public:
    int off[3];
    OffCount (void) {
        off[0] = 0;
        off[1] = 0;
        off[2] = 0;
    }
    OffCount (int num) {
        off[0] = 0;
        off[1] = 0;
        off[2] = 0;
        off[num] = 1;
    }
};

map< string, vector<Info_Transcript*>* > transc_per_chrom_dic;
map<string, int> seq_pk_dic;
map<string, OffCount> cnt_dic;

mutex m_fo;
mutex m_cnt_dic;
mutex m_nf;
mutex m_file_cnt;
semaphore sem_thread(MAX_THREAD_CNT);

int file_cnt = 0;

int pk_off = 1;

void run_thread(string &p, ofstream &info_file, ofstream &fo1) {
    int num_off, transcid;
    ostringstream ss;

    sem_thread.wait();

    map<string, int>::iterator it2;
    map<string, OffCount>::iterator it3;

    vector<string> entries;
    string line;

    m_file_cnt.lock();
    ss << p << "/outs_" << file_cnt+1 << ".txt";
    file_cnt++;
    m_file_cnt.unlock();

    ifstream f(ss.str().c_str());
    while (getline(f, line)) {
        entries = split(rstrip(line), '\t');
        num_off = stoi(entries[5]);
        entries[0].resize(20);
        it2 = seq_pk_dic.find(entries[0]);
        if (it2 == seq_pk_dic.end()) { // Not found
            m_nf.lock();
            fo1 << entries[0] << endl;
            m_nf.unlock();
        } else {
            if (transc_per_chrom_dic.find(entries[1]) == transc_per_chrom_dic.end()) continue;
            transcid = find_transcript(transc_per_chrom_dic[entries[1]], stoi(entries[2]), entries[4][0]);
            m_fo.lock();
            info_file << pk_off << '\t'
                      << (*it2).second << '\t'
                      << num_off << '\t'
                      << entries[1] << '\t'
                      << entries[2] << '\t'
                      << entries[3] << '\t'
                      << entries[4] << '\t'
                      << transcid << endl;
            m_cnt_dic.lock();
            it3 = cnt_dic.find(entries[0]);
            if (it3 == cnt_dic.end()) { // Not found
                OffCount *cnt = new OffCount(num_off);
                cnt_dic[entries[0]] = *cnt;
            } else {
                cnt_dic[entries[0]].off[num_off]++;
            }
            m_cnt_dic.unlock();
            pk_off++;
            m_fo.unlock();
        }
    }
    f.close();
    sem_thread.notify();
}

int main(int argc, char** argv) {
    vector<string> entries;
    string p = argv[1];
    p = rstrip(p, '/');
    
    cout << "Reading required files..." << endl;
    string line;
    ifstream f1((p+"/info_guides.txt").c_str());
    while (getline(f1, line)) {
        entries = split(rstrip(line), '\t');
        seq_pk_dic[entries[1]] = stoi(entries[0]);
    }
    f1.close();

    ifstream f2((p+"/info_transcripts.txt").c_str());
    vector<Info_Transcript*> *transcripts;
    map< string, vector<Info_Transcript*>* >::iterator it1;
    while (getline(f2, line)) {
        entries = split(rstrip(line), '\t');
        it1 = transc_per_chrom_dic.find(entries[2]);
        if (it1 == transc_per_chrom_dic.end()) {  // Not found
            transc_per_chrom_dic[entries[2]] = new vector<Info_Transcript*>;
            transcripts = transc_per_chrom_dic[entries[2]];
        } else {
            transcripts = (*it1).second;
        } 
        transcripts->push_back(new Info_Transcript(stoi(entries[0]), slice(entries, 4)));
    }
    f2.close();

    ofstream info_file((p+"/info_cpp_offtargets_ref.txt").c_str());

    cout << "Counting off-targets and locating its origin..." << endl;
    int i;
    ofstream fo1((p+"/offtarget_cpp_notfound.txt").c_str());
    thread *threads[FILE_CNT];
    for (i=0; i<FILE_CNT; i++) {
        threads[i] = new thread([&] { run_thread(p, info_file, fo1); });
    }
    for (i=0; i<FILE_CNT; i++) threads[i]->join();
    info_file.close();
    fo1.close();

    cout << "Printing results..." << endl;
    ofstream fo2((p+"/info_cpp_offtargets.txt").c_str());
    typedef map<string, OffCount>::iterator it_type;
    for (it_type iterator = cnt_dic.begin(); iterator != cnt_dic.end(); iterator++) {
        fo2 << seq_pk_dic[iterator->first] << '\t'
            << iterator->second.off[0] << '\t'
            << iterator->second.off[1] << '\t'
            << iterator->second.off[2] << endl;
    }
    return 0;
}
