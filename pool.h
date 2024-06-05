/*************************************************************************
    > File Name: pool.h
    > Author: Haoming Bai
    > Mail: haomingbai@hotmail.com
    > Created Time: 2024年05月23日 星期四 22时46分04秒
 ************************************************************************/

class memory_block {
private:
	char *start, *next, *prev;
	int count;
public:
	memory_block():count(0) {
		start = new char[4096];
	}
	~memory_block() {
		delete[] (char *)start;
	}
};

class memory_pool {};
