/*
 * IFlowCutter.cpp
 *   This file is a degraded copy of pace.cpp in flow-cutter-pace17.
 *   I made a small change to the source code for our use (e.g., timeout).
 *
 *  Editted on: 2022/03/29
 *      Author: k-hasimt
 */


//#include "flow-cutter-pace17/src/id_func.h"
//#include "flow-cutter-pace17/src/list_graph.h"
#include "flow-cutter-pace17/src/multi_arc.h"
//#include "flow-cutter-pace17/src/sort_arc.h"
//#include "flow-cutter-pace17/src/chain.h"
#include "flow-cutter-pace17/src/union_find.h"
//#include "flow-cutter-pace17/src/node_flow_cutter.h"
#include "flow-cutter-pace17/src/separator.h"
#include "flow-cutter-pace17/src/id_multi_func.h"
#include "flow-cutter-pace17/src/filter.h"
#include "flow-cutter-pace17/src/preorder.h"
#include "flow-cutter-pace17/src/contraction_graph.h"
//#include "flow-cutter-pace17/src/tree_decomposition.h"
#include "flow-cutter-pace17/src/greedy_order.h"
//#include "flow-cutter-pace17/src/min_max.h"
#include "flow-cutter-pace17/src/heap.h"

#include "preprocessor/IFlowCutter.h"
#include "preprocessor/TreeDecomposition.h"
#include "utils/System.h"

// #include <limits>
// #include <signal.h>
// #include <stdlib.h>
// #include <string.h>
// #include <string>
// #include <sstream>
#include <queue>

#include <sys/time.h>
// #include <unistd.h>

#include <iostream>
using namespace std;


IFlowCutter::IFlowCutter(int n, int m, double timeout) :
		nodes(n), best_bag_size(numeric_limits<int>::max()), head(2*m, n), tail(2*m, n), timeout(timeout) {}

void IFlowCutter::importGraph(const Graph& g)
{
	int next_arc = 0;
	for(int i=0; i<nodes; i++)
		for(auto j : g.Neighbors(i)) {
			if(i >= j) continue;
			assert(next_arc < head.preimage_count());
			head[next_arc] = j;
			tail[next_arc] = i;
			next_arc++;
			head[next_arc] = i;
			tail[next_arc] = j;
			next_arc++;
		}
}

void ignore_return_value(long long){}

unsigned long long get_milli_time(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (unsigned long long)(tv.tv_sec) * 1000
	   + (unsigned long long)(tv.tv_usec) / 1000;
}

// This hack is actually standard compilant
template <class T, class S, class C>
S& access_internal_vector(std::priority_queue<T, S, C>& q) {
	struct Hacked : private priority_queue<T, S, C> {
		static S& access(priority_queue<T, S, C>& q) {
			return q.*&Hacked::c;
		}
	};
	return Hacked::access(q);
}

void print_comment(std::string msg){
	msg = "c o "+std::move(msg) + "\n";
	ignore_return_value(write(STDOUT_FILENO, msg.data(), msg.length()));
}

template<class Tail, class Head>
void check_multilevel_partition_invariants(const Tail&tail, const Head&head, const std::vector<Cell>&multilevel_partition){
	#ifndef NDEBUG
	const int node_count = tail.image_count();
	const int arc_count = tail.preimage_count();

	auto is_child_of = [&](int c, int p){
		for(;;){
			if(c == p)
				return true;
			if(c == -1)
				return false;
			c = multilevel_partition[c].parent_cell;
		}
	};

	auto are_ordered = [&](int a, int b){
		return is_child_of(a, b) || is_child_of(b, a);
	};

	ArrayIDFunc<int> cell_of_node(node_count);
	cell_of_node.fill(-1);

	for(int i=0; i<(int)multilevel_partition.size(); ++i){
		for(auto&y:multilevel_partition[i].separator_node_list){
			assert(cell_of_node(y) == -1);
			cell_of_node[y] = i;
		}
	}

	for(auto x:cell_of_node)
		assert(x != -1);

	for(int xy = 0; xy < arc_count; ++xy){
		int x = cell_of_node(tail(xy)), y = cell_of_node(head(xy));
		assert(are_ordered(x, y));
	}
	#endif
}

template<class Tail, class Head, class ComputeSeparator, class OnNewMP>
void compute_multilevel_partition(const Tail&tail, const Head&head, const ComputeSeparator&compute_separator, int smallest_known_treewidth, const OnNewMP&on_new_multilevel_partition){

	const int node_count = tail.image_count();
	const int arc_count = tail.preimage_count();

	std::vector<Cell>closed_cells;
	std::priority_queue<Cell>open_cells;

	{
		Cell top_level_cell;
		top_level_cell.separator_node_list.resize(node_count);
		for(int i=0; i<node_count; ++i)
			top_level_cell.separator_node_list[i] = i;
		//top_level_cell.boundary_node_list = {};
		top_level_cell.parent_cell = -1;

		open_cells.push(std::move(top_level_cell));
	}

	int max_closed_bag_size = 0;
	int max_open_bag_size = node_count;

	auto check_if_better = [&]{
		int current_tree_width = std::max(max_closed_bag_size, max_open_bag_size);

		if(current_tree_width < smallest_known_treewidth){
			smallest_known_treewidth = current_tree_width;

			std::vector<Cell>cells = closed_cells;
			for(auto&q:access_internal_vector(open_cells))
				cells.push_back(q);
			check_multilevel_partition_invariants(tail, head, cells);
			on_new_multilevel_partition(cells, open_cells.empty() || max_closed_bag_size>=max_open_bag_size);
		}
	};

	check_if_better();


	ArrayIDFunc<int>node_to_sub_node(node_count);
	node_to_sub_node.fill(-1);

	auto inv_tail = invert_sorted_id_id_func(tail);

	BitIDFunc in_child_cell(node_count);
	in_child_cell.fill(false);

	while(!open_cells.empty()){

		#ifndef NDEBUG

		int real_max_closed_bag_size = 0;
		for(auto&x:closed_cells)
			max_to(real_max_closed_bag_size, x.bag_size());
		assert(max_closed_bag_size == real_max_closed_bag_size);

		int real_max_open_bag_size = 0;
		for(auto&x:access_internal_vector(open_cells))
			max_to(real_max_open_bag_size, x.bag_size());
		assert(max_open_bag_size == real_max_open_bag_size);

		#endif

		auto current_cell = std::move(open_cells.top());
		open_cells.pop();

		bool must_recompute_max_open_bag_size = (current_cell.bag_size() == max_open_bag_size);

		int closed_cell_id = closed_cells.size();

		if(current_cell.bag_size() > max_closed_bag_size){

			auto interior_node_list = std::move(current_cell.separator_node_list);
			int interior_node_count = interior_node_list.size();

			ArrayIDFunc<int>sub_node_to_node(interior_node_count);

			int next_sub_id = 0;
			for(int x:interior_node_list){
				node_to_sub_node[x] = next_sub_id;
				sub_node_to_node[next_sub_id] = x;
				++next_sub_id;
			}

			auto is_node_interior = id_func(
				node_count,
				[&](int x)->bool{
					return node_to_sub_node(x) != -1;
				}
			);

			auto is_arc_interior = id_func(
				arc_count,
				[&](int xy)->bool{
					return is_node_interior(tail(xy)) && is_node_interior(head(xy));
				}
			);

			int interior_arc_count = count_true(is_arc_interior);
			auto sub_tail = keep_if(is_arc_interior, interior_arc_count, tail);
			auto sub_head = keep_if(is_arc_interior, interior_arc_count, head);

			for(auto&x:sub_tail)
				x = node_to_sub_node(x);
			sub_tail.set_image_count(interior_node_count);

			for(auto&x:sub_head)
				x = node_to_sub_node(x);
			sub_head.set_image_count(interior_node_count);

			auto sub_separator = compute_separator(sub_tail, sub_head);

			BitIDFunc is_in_sub_separator(interior_node_count);
			is_in_sub_separator.fill(false);
			for(auto x:sub_separator)
				is_in_sub_separator.set(x, true);

			UnionFind uf(interior_node_count);

			for(int xy=0; xy<interior_arc_count; ++xy){
				int x = sub_tail(xy);
				int y = sub_head(xy);
				if(!is_in_sub_separator(x) && !is_in_sub_separator(y))
					uf.unite(x, y);
			}

			std::vector<std::vector<int>>nodes_of_representative(interior_node_count);
			for(int x=0; x<interior_node_count; ++x)
				if(!is_in_sub_separator(x))
					nodes_of_representative[uf(x)].push_back(x);

			auto&separator = sub_separator;
			for(auto&x:separator)
				x = sub_node_to_node(x);

			for(int x=0; x<interior_node_count; ++x){
				if(!nodes_of_representative[x].empty()){
					Cell new_cell;

					auto&new_cell_interior_node_list = nodes_of_representative[x];
					for(auto&x:new_cell_interior_node_list)
						x = sub_node_to_node(x);

					new_cell.parent_cell = closed_cell_id;

					new_cell.separator_node_list = std::move(new_cell_interior_node_list);

					new_cell.boundary_node_list = current_cell.boundary_node_list;
					new_cell.boundary_node_list.insert(new_cell.boundary_node_list.end(), separator.begin(), separator.end());

					{
						for(auto x:new_cell.separator_node_list)
							in_child_cell.set(x, true);
						new_cell.boundary_node_list.erase(
							std::remove_if(
								new_cell.boundary_node_list.begin(),
								new_cell.boundary_node_list.end(),
								[&](int x)->bool{
									for(auto xy:inv_tail(x))
										if(in_child_cell(head(xy)))
											return false;
									return true;
								}
							),
							new_cell.boundary_node_list.end()
						);
						for(auto x:new_cell.separator_node_list)
							in_child_cell.set(x, false);
					}

					new_cell.separator_node_list.shrink_to_fit();
					new_cell.boundary_node_list.shrink_to_fit();

					if(new_cell.bag_size() > max_open_bag_size)
						max_open_bag_size = new_cell.bag_size();

					open_cells.push(std::move(new_cell));
				}
			}

			current_cell.separator_node_list = std::move(separator);
			current_cell.separator_node_list.shrink_to_fit();

			for(int x:interior_node_list)
				node_to_sub_node[x] = -1;
		}

		if(current_cell.bag_size() > max_closed_bag_size)
			max_closed_bag_size = current_cell.bag_size();

		if(must_recompute_max_open_bag_size){
			max_open_bag_size = 0;
			for(auto&x:access_internal_vector(open_cells))
				if(x.bag_size() > max_open_bag_size)
					max_open_bag_size = x.bag_size();
		}

		closed_cells.push_back(std::move(current_cell));

		check_if_better();

		if(max_closed_bag_size >= smallest_known_treewidth){
			return;
		}

		if(max_closed_bag_size >= max_open_bag_size){
			return;
		}
	}
}

int IFlowCutter::compute_max_bag_size_of_order(const ArrayIDIDFunc&order){
	auto inv_order = inverse_permutation(order);
	int current_tail = -1;
	int current_tail_up_deg = 0;
	int max_up_deg = 0;
	compute_chordal_supergraph(
		chain(tail, inv_order), chain(head, inv_order),
		[&](int x, int y){
			if(current_tail != x){
				current_tail = x;
				max_to(max_up_deg, current_tail_up_deg);
				current_tail_up_deg = 0;
			}
			++current_tail_up_deg;
		}
	);
	return max_up_deg+1;
}

void IFlowCutter::test_new_order(const ArrayIDIDFunc&order, TreeDecomposition& td){
	int x = compute_max_bag_size_of_order(order);
	{
		if(x < best_bag_size){
			best_bag_size = x;
			td = output_tree_decompostion_of_order(tail, head, order);
		}
	}
}

TreeDecomposition IFlowCutter::output_tree_decompostion_of_order(ArrayIDIDFunc tail, ArrayIDIDFunc head, const ArrayIDIDFunc&order){
	TreeDecomposition better_td;

	const int node_count = tail.image_count();

	auto inv_order = inverse_permutation(order);
	tail = chain(tail, inv_order);
	head = chain(head, inv_order);

	vector<vector<int>>nodes_in_bag;
	ArrayIDFunc<vector<int>>bags_of_node(node_count);

	auto is_left_subset_of_right = [](const std::vector<int>&l, const std::vector<int>&r){
		auto i = l.begin(), j = r.begin();

		for(;;){
			if(i == l.end())
				return true;
			if(j == r.end())
				return false;

			if(*i < *j)
				return false;
			if(*i == *j)
				++i;
			++j;
		}
	};

	auto compute_intersection_size = [](const std::vector<int>&l, const std::vector<int>&r){
		auto i = l.begin(), j = r.begin();
		int n = 0;
		for(;;){
			if(i == l.end() || j == r.end())
				return n;
			if(*i < *j)
				++i;
			else if(*i > *j)
				++j;
			else{
				++i;
				++j;
				++n;
			}
		}
	};

	auto on_new_potential_maximal_clique = [&](int lowest_node_in_clique, std::vector<int>clique){
		for(auto b:bags_of_node(lowest_node_in_clique))
			if(is_left_subset_of_right(clique, nodes_in_bag[b]))
				return;
		int bag_id = nodes_in_bag.size();
		for(auto x:clique)
			bags_of_node[x].push_back(bag_id);
		nodes_in_bag.push_back(move(clique));
	};


	{
		BitIDFunc is_root(node_count);
		is_root.fill(true);
		std::vector<int>upper_neighborhood_of_z;
		int z = -1;
		compute_chordal_supergraph(
			tail, head,
			[&](int x, int y){
				is_root.set(x, false);
				if(z != -1 && z != x){
					upper_neighborhood_of_z.push_back(z);
					sort(upper_neighborhood_of_z.begin(), upper_neighborhood_of_z.end());
					on_new_potential_maximal_clique(z, move(upper_neighborhood_of_z));
					upper_neighborhood_of_z.clear();
				}
				z = x;
				upper_neighborhood_of_z.push_back(y);
			}
		);
		if(z != -1){
			upper_neighborhood_of_z.push_back(z);
			sort(upper_neighborhood_of_z.begin(), upper_neighborhood_of_z.end());
			on_new_potential_maximal_clique(z, move(upper_neighborhood_of_z));
		}

		for(int x=0; x<node_count; ++x){
			if(is_root(x)){
				on_new_potential_maximal_clique(x, {x});
			}
		}
	}



	int bag_count = nodes_in_bag.size();

	int maximum_bag_size = 0;
	for(auto&b:nodes_in_bag)
		if((int)b.size() > maximum_bag_size)
			maximum_bag_size = b.size();

	better_td.init(bag_count);
	better_td.setWidth(maximum_bag_size-1);
	better_td.setNumGraphNodes(node_count);
	cout << "c o td: #bags " << bag_count << ", tw " << maximum_bag_size
			<< ", elapsed " << Glucose::cpuTime()  << " s"<< endl; // << endl;// << ", #vars " << node_count << endl;
	better_td.initBags();

	for(int i=0; i<bag_count; ++i) {
		vector<vector<int>>& bags = better_td.Bags();
		for(auto x: nodes_in_bag[i])
			bags[i].push_back(order(x));
	}

	{
		std::vector<int>tail, head, weight;

		for(int b=0; b<bag_count; ++b){
			vector<int>neighbor_bags;
			for(auto x:nodes_in_bag[b]){
				vector<int>tmp;
				std::set_union(
					bags_of_node[x].begin(), bags_of_node[x].end(),
					neighbor_bags.begin(), neighbor_bags.end(),
					std::back_inserter(tmp)
				);
				neighbor_bags.swap(tmp);
			}
			for(auto p:neighbor_bags){
				if(p != b){
					tail.push_back(b);
					head.push_back(p);
					weight.push_back(compute_intersection_size(nodes_in_bag[b], nodes_in_bag[p]));
				}
			}
		}

		int arc_count = tail.size();

		auto out_arc = invert_id_id_func(
			id_id_func(
				arc_count, bag_count,
				[&](unsigned a){return tail[a];}
			)
		);

		BitIDFunc in_tree(bag_count);
		in_tree.fill(false);
		max_id_heap<int>q(arc_count);

		for(int b=0; b<bag_count; ++b){
			if(!in_tree(b)){
				if(b != 0)
					better_td.addEdge(0, b);
				in_tree.set(b, true);
				for(int a:out_arc(b))
					q.push(a, weight[a]);
				while(!q.empty()){
					int xy = q.pop();

					int x = tail[xy];
					int y = head[xy];

					assert(in_tree(x));

					if(!in_tree(y)){
						better_td.addEdge(x, y);
						in_tree.set(y, true);
						for(int yz:out_arc(y)){
							assert(!q.contains(yz));
							q.push(yz, weight[yz]);
						}
					}
				}

			}
		}
	}

	return better_td;
}

TreeDecomposition IFlowCutter::output_tree_decompostion_of_multilevel_partition(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head, const ArrayIDIDFunc&to_input_node_id, const std::vector<Cell>&cell_list){
	TreeDecomposition better_td;

	int bag_count = cell_list.size();

	better_td.init(bag_count);
	better_td.setWidth(get_treewidth_of_multilevel_partition(cell_list)-1);
	better_td.setNumGraphNodes(get_node_count_of_multilevel_partition(cell_list));
	cout << "c o td: #bags " << bag_count << ", tw " << get_treewidth_of_multilevel_partition(cell_list)-1
			<< ", elapsed " << Glucose::cpuTime()  << " s"<< endl; // << ", #vars " << get_node_count_of_multilevel_partition(cell_list) << endl;
	better_td.initBags();

	for(int i=0; i<bag_count; ++i) {
		vector<vector<int>>& bags = better_td.Bags();
		for(auto&x:cell_list[i].separator_node_list)
			bags[i].push_back(to_input_node_id(x));
		for(auto&x:cell_list[i].boundary_node_list)
			bags[i].push_back(to_input_node_id(x));
	}

	for(int i=0; i<bag_count; ++i){
		if(cell_list[i].parent_cell != -1)
			better_td.addEdge(i, cell_list[i].parent_cell);
	}

	return better_td;
}

TreeDecomposition IFlowCutter::constructTD()
{
	double start_time = Glucose::cpuTime();

	TreeDecomposition td;

	ArrayIDIDFunc preorder, inv_preorder;

	int random_seed = 0;
	try{
		{
			preorder = compute_preorder(compute_successor_function(tail, head));
			for(int i=0; i<tail.image_count(); ++i)
				preorder[i] = i;
			inv_preorder = inverse_permutation(preorder);
			tail = chain(std::move(tail), inv_preorder);
			head = chain(std::move(head), inv_preorder);
		}

		{
			auto p = sort_arcs_first_by_tail_second_by_head(tail, head);
			tail = chain(p, std::move(tail));
			head = chain(p, std::move(head));
		}

		const int node_count = tail.image_count();

		long long last_print = 0;

		auto on_new_multilevel_partition = [&](const std::vector<Cell>&multilevel_partition, bool must_print){

			long long now = get_milli_time();

			if(!must_print && now - last_print < 30000)
				return;
			last_print = now;

			int tw = get_treewidth_of_multilevel_partition(multilevel_partition);
			{
				/* update best tree decomposition*/

				td = output_tree_decompostion_of_multilevel_partition(tail, head, preorder, multilevel_partition);
				best_bag_size = tw;
			}
			// print_comment("status "+to_string(best_bag_size)+" "+to_string(get_milli_time()));
			// print_comment("status "+to_string(best_bag_size)+" "+to_string(get_milli_time()));
			start_time = Glucose::cpuTime();
		};


		{
			try{
				std::minstd_rand rand_gen;
				rand_gen.seed(random_seed);

				if(node_count > 500000)
				{
					print_comment("start F1 with 0.1 min balance and edge_first");
					flow_cutter::Config config;
					config.cutter_count = 1;
					config.random_seed = rand_gen();
					config.min_small_side_size = 0.1;
					config.max_cut_size = 500;
					config.separator_selection = flow_cutter::Config::SeparatorSelection::edge_first;
					compute_multilevel_partition(tail, head, flow_cutter::ComputeSeparator(config), best_bag_size, on_new_multilevel_partition);
				}

				if(node_count < 50000){
					print_comment("min degree heuristic");
					test_new_order(chain(compute_greedy_min_degree_order(tail, head), inv_preorder), td);
				}

				if(node_count < 10000){
					print_comment("min shortcut heuristic");
					test_new_order(chain(compute_greedy_min_shortcut_order(tail, head), inv_preorder), td);
				}

				start_time = Glucose::cpuTime();
				{
					print_comment("run with 0.0/0.1/0.2 min balance and node_min_expansion in endless loop with varying seed");
					flow_cutter::Config config;
					config.cutter_count = 1;
					config.random_seed = rand_gen();
					config.max_cut_size = 10000;
					config.separator_selection = flow_cutter::Config::SeparatorSelection::node_min_expansion;

					for(int i=2;;++i){
						config.random_seed = rand_gen();
						if(i % 16 == 0)
							++config.cutter_count;

						switch(i % 3){
						case 2: config.min_small_side_size = 0.2; break;
						case 1: config.min_small_side_size = 0.1; break;
						case 0: config.min_small_side_size = 0.0; break;
						}

						compute_multilevel_partition(tail, head, flow_cutter::ComputeSeparator(config), best_bag_size, on_new_multilevel_partition);

						if(Glucose::cpuTime() - start_time > timeout) break;
					}
				}
			}catch(...){
			}

		}
	}catch(...){
	}

	return td;
}

