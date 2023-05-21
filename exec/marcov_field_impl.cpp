// marcov_field_impl.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
#include "quad_tree_handler.h"
#include <boost/filesystem.hpp>
#include "optional_handler.h"
#include <fstream>
#include "network_prob_img.h"
#include "network_prob_img_ensemble.h"
#include <json/json.h>

using namespace boost::filesystem;
//int main(){
int main(int argc, char *argv[]) {

	//тест для обработки большого реального изображения после прохода нейросетью
	optional_handler t;

	// дефолтные значения аргументов для запуска обработчика
	//"C:\\Users\\anastasya\\Desktop\\data_FR.txt", false);
	//("C:\\Users\\anastasya\\Desktop\\data_test.txt", false);
	string fileDescrName = "C:\\Users\\anastasya\\Desktop\\data_Kubinka.txt", fileResName = "D:\\_SAR_Kubinka\\class_results.txt";
	string fileFaultClassName = "D:\\_SAR_Kubinka\\fault_results.txt", fileClassMarks = "D:/_SAR_Kubinka/class_marks.txt";
	string modelName = "", csvFileName = "";
	basic_prob_img *img;
	int handlerType = 5, handler_mode = 2, img_mode = 2, class_amount = 4, layer_order = 2;
	float cut_pres = 0.0;
	vector<string> fileNames/* = {"D:/_SAR_Kubinka/class_probs_0.txt",
"D:/_SAR_Kubinka/class_probs_1.txt",
"D:/_SAR_Kubinka/class_probs_2.txt",
"D:/_SAR_Kubinka/class_probs_3.txt" }*/;
	string borders_filename = "";
	cout << "in programm" << endl;
	if (argc != 1)
	{
		std::ifstream load_params;
		int files_amount;
		string buf;
		
		Json::Value root;
		Json::Reader reader;
		reader.parse(argv[1], root);
		
		handlerType  = root["worker_number"].asInt();
		handler_mode = root["handler_mode"].asInt();
		img_mode     = root["img_mode_gener"].asInt();
		class_amount = root["class_amount"].asInt();
		modelName    = root["model"].asString();
		layer_order  = root["layer_ordering"].asInt();
		cut_pres     = root["cut_presision"].asFloat();
		csvFileName  = root["csv_filename"].asString();
		
		load_params.open(root["config_file"].asString());
		// считываем названия основных конфигурационных файлов для обработки результатов нейросети
		load_params >> files_amount;
		load_params >> buf;
		fileDescrName = string(buf.begin(), buf.end());
		cout << fileDescrName << endl;
		load_params >> buf;
		fileClassMarks = string(buf.begin(), buf.end());
		load_params >> buf;
		fileResName = string(buf.begin(), buf.end());
		load_params >> buf;
		fileFaultClassName = string(buf.begin(), buf.end());
		// здесь считываем информацию о результатах классифкации обычных и накопленных изображений,  сделанных нейросеткой
		load_params >> files_amount;
		for (int i = 0; i < files_amount; ++i)
		{
			load_params >> buf;
			fileNames.push_back(string(buf.begin(), buf.end()));
		}
		load_params >> borders_filename;
		load_params.close();
	}
	
	switch (handlerType) 
	{
	case 1:
	{
		img = new initial_prob_img();
		img->gen_prob_img_from_config(fileDescrName, img_mode);
		shared_ptr<basic_prob_img> ptr_img(img);
		cout << "into mixture" << endl;
		t.mixture_handler(ptr_img->get_m_image(), class_amount, handler_mode, 0.001);
	}
	break;
	case 2:
	{
		img = new initial_prob_img();  
		img->gen_prob_img_from_config(fileDescrName, img_mode);
		shared_ptr<basic_prob_img> ptr_img(img);
		t.quadtree_handler(ptr_img, class_amount, cut_pres, layer_order, handler_mode);
	}
	break;
	case 3:
	{
		
		img = new newtwork_prob_img();
		
		img->gen_prob_img_from_config(fileDescrName, img_mode);
		
		shared_ptr<basic_prob_img> ptr_img(img);
		
		t.network_results_handler(ptr_img->get_m_image(), class_amount, handler_mode, fileClassMarks);
	}
	break;
	case 4:
	{
		
		img = new newtwork_prob_img();
		img->gen_prob_img_from_config(fileDescrName, img_mode);
		img->load_probs_from_file(fileNames);
		shared_ptr<basic_prob_img> ptr_img(img);	
		t.quadtree_handler(ptr_img, class_amount, cut_pres, layer_order, handler_mode);
		
	}
	break;
	case 5:
	{

		img = new newtwork_prob_img_ensemble();
		img->gen_prob_img_from_config(fileDescrName, img_mode);
		cout << "gen " << fileNames.size() << endl;
		img->load_probs_from_file(fileNames, borders_filename);
		shared_ptr<basic_prob_img> ptr_img(img);
		//t.quadtree_handler(ptr_img, 5, 0.0, 1);
		//cout << "gen2" << endl;
		t.quadtree_handler_ensemble(ptr_img, class_amount, cut_pres, layer_order, handler_mode);
		//t.network_results_handler(ptr_img->get_m_image(), 5, fileClassMarks);
	}
	break;
	case 6:
	{

		img = new newtwork_prob_img_ensemble();
		img->gen_prob_img_from_config(fileDescrName, img_mode);
		cout << "gen " << fileNames.size() << endl;
		img->load_probs_from_file(fileNames, borders_filename);
		shared_ptr<basic_prob_img> ptr_img(img);
		
		t.quadtree_handler_ensemble_sqip(ptr_img, class_amount, cut_pres, layer_order, handler_mode);
		//t.network_results_handler(ptr_img->get_m_image(), 5, fileClassMarks);
	}
	break;
	case 7:
	{
		img = new initial_prob_img();
		img->gen_prob_img_from_config(fileDescrName, img_mode);
		shared_ptr<basic_prob_img> ptr_img(img);
		t.extended_image_handler(ptr_img, class_amount, cut_pres, layer_order, handler_mode);
	}
	break;
	}


	//t.detect_result_by_mask(fileResName, fileFaultClassName, modelName, csvFileName);
	t.printInformation();
	//t.draw_graphics();



	////тест для обработки маленького собственно-сгенеренного изображения
	//{
	//	initial_prob_img *img = new initial_prob_img(32, NORMAL, 1, 2);
	////initial_prob_img *img = new initial_prob_img(32, RAYLEIGH, 1, 2);
	//	shared_ptr<initial_prob_img> ptr_img(img);
	//	quad_tree_handler bbbb = quad_tree_handler(ptr_img, 32, nullptr);
	//	bbbb.set_probabilities(0, 0);
	//	bbbb.bottom_up_pass();
	//	auto begin1 = std::chrono::steady_clock::now();
	//	bbbb.up_down_pass_V2();
	//	//bbbb.up_down_pass();
	//	auto end1 = std::chrono::steady_clock::now();
	//	auto elapsed_ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1);
	//	cout << "elapsed_ms1  " << elapsed_ms1.count() << endl;
	//	////bbbb.split_image_by_summ();
	//	bbbb.split_image_by_vote();
	//	bbbb.create_splitted_img();
	//	bbbb.draw_graphics();
	//}
	////


	////тест для обработки маленького собственно-сгенеренного изображения  по нескольку раз одного и того же 
	//{
	//	/*initial_prob_img *img = new initial_prob_img("D:\\generated_image.txt", 32, RAYLEIGH, 1, 2);*/
	//	//initial_prob_img *img = new initial_prob_img("D:\\generated_image.txt", 32, NORMAL, 1, 2);
	//	shared_ptr<initial_prob_img> ptr_img(img);
	//	quad_tree_handler bbbb = quad_tree_handler(ptr_img, 32);
	//	bbbb.set_probabilities(0, 0);
	//	bbbb.bottom_up_pass();
	//	bbbb.up_down_pass();
	//	////bbbb.split_image_by_summ();
	//	bbbb.split_image_by_vote();
	//	bbbb.create_splitted_img();
	//	bbbb.draw_graphics();
	//}
	////
	system("pause");
	return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
