#pragma once
#include "mix_img_obj.h"
#include "basic_prob_img.h"
#include "omp.h"

class newtwork_prob_img: public basic_prob_img
{

public:
    newtwork_prob_img(){}

    // ����������� �������� ��� ����� � �����������, � ����� ����, ����� ��� ��� �������������� �������� �����������
    virtual void gen_prob_img_from_config(string filename, int mode);//, bool flag);
    virtual void load_probs_from_file(vector<string> probs_data, string borders_file = "");

};


