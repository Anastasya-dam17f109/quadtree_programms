#pragma once
#include "mix_img_obj.h"
#include "omp.h"
#include <vector>
#include <unordered_map>
#include <cstdlib>

using namespace std;


struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const
    {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);

        if (hash1 != hash2) {
            return hash1 ^ hash2;
        }

        // If hash1 == hash2, their XOR is zero.
        return hash1;
    }
};

class basic_prob_img
{
public:
    shared_ptr<mix_img_obj> m_image;
	// ������� ����� �����������
    shared_ptr <int[]>      m_layer_size;
    shared_ptr <int[]>      m_layer_idx;
	// ������� ��������� ������ ����� � ���������� ������� ������������
    shared_ptr <int[]>      m_init_layer_idx;
    std::unordered_map<std::pair<int, int>, int, hash_pair> m_homo_areas;
    int m_layer_amount = 1;
	// ���������� ������ ������������
    double* m_prob_img;

    shared_ptr<double[]> mix_shift;
    shared_ptr<double[]> mix_scale;
    shared_ptr<double[]> m_cl_probs;
    unsigned m_image_len_x = 32;
    unsigned m_image_len_y = 32;
    unsigned m_class_amount = 1;
    int m_init_window_size = 8;

    basic_prob_img(){}
	// ��������� ������� ����������� � ������������� ������������ ������ (� ��� ����������)
    virtual void gen_prob_img_with_targs(int img_size, mix_type mix_t, int amount_targets, int classes){}
	// ��������� ������� ����������� � ������������� ������������ ������ (� ��� ��������� �� �����)
    virtual void gen_prob_img_with_targs_from_file(string file_name, int img_size, mix_type mix_t, int amount_targets, int classes){}
	// ��������� ������� ����������� , ����������� � ������ 
    virtual void gen_prob_img_from_obj(shared_ptr<mix_img_obj> image){}

	// ���������� ������������ ����� �����,  ���������� � ����� ������������
    virtual void gen_prob_img_from_config(string filename, int mode){}
    virtual void load_probs_from_file(vector<string> probs_data, string borders_file="") {}
	// ��������� ������ ��� ������ ������������
    void  alloc_layer_mmr()
    {
        m_init_layer_idx = shared_ptr <int[]>(new int [m_layer_amount]);
        int summ = 0;
		
        for (int k = 0; k < m_layer_amount; ++k) {
            m_init_layer_idx[k] = summ;
			
            summ += m_layer_size[2*k] * m_layer_size[2*k+1] * m_class_amount;
        }
        m_prob_img = new double [summ];
    }

    int                     get_class_amount()    {return m_class_amount;}
    std::pair<int, int>     get_image_len()       {return std::pair<int, int>(m_image_len_x, m_image_len_y); }
    double*                 get_image()           {return m_prob_img;}
    shared_ptr<mix_img_obj> get_m_image()         {return m_image;}
    shared_ptr<double[]>    get_cl_probs()        {return m_cl_probs;}
    int                     get_layer_amount()    {return m_layer_amount;}
    shared_ptr <int[]>      get_init_layer_idx () {return m_init_layer_idx;}
    shared_ptr <int[]>      get_init_layer_size() {return m_layer_size;}
    std::unordered_map<std::pair<int, int>, int, hash_pair>& get_homo_areas() { return m_homo_areas; }
    virtual ~basic_prob_img()
    {
        delete [] m_prob_img;
    }
};

