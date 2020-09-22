#include "TSP.h"

using namespace std;

bool compareL1(Individuo i1, Individuo i2){
	return (i1.L1() < i2.L1()); 
} 

void City::set_pos(int & posiz){
	pos = posiz;
}

Individuo::Individuo(int number, const City & inizio, vector<City>::iterator it1, vector<City>::iterator it2) : N(number), partenza(inizio) {
	for (auto it = it1; it !=it2; it++) {
		percorso.push_back(*it);
	}
}

//Scelgo casualmente due delle ultime 31 città (da 0 a 30) e le scambio. Pesco da 0 a 31 così la parte intera sarà da 0 a 30 e sono 31 valori.
void Individuo::mutation_swap(Random & rnd){
	int i = int(rnd.Rannyu(0,31));
	int j=0;
	do{
       		j = int(rnd.Rannyu(0,31));
     	}while(i==j);
	swap(this->percorso[i],this->percorso[j]);
	rnd.SaveSeed();
}

//Scelgo casualmente un numero n<N-1 e traslo tutte le città (tranne la prima) di n.
void Individuo::mutation_shift(Random & rnd){
	int n = int(rnd.Rannyu(0,30));
	vector<City> copia; 
	for (auto it = this->percorso.begin()+int(this->percorso.size()-n); it != this->percorso.end(); it++) {
		copia.push_back(*it);
	}
	for (auto it = this->percorso.begin(); it != this->percorso.begin()+(percorso.size()-n); it++) {
		copia.push_back(*it);
	}
	int m  = this->percorso.size();
	for (int i=0; i< m; i++){
		swap(this->percorso[i],copia[i]);
	}	
	copia.clear();
	rnd.SaveSeed();
}

//Scelgo casualmente un numero m<N/2 e permuto m città contigue con altre m città contigue.
void Individuo::mutation_permutation(Random & rnd){
	int m = int(rnd.Rannyu(0,15));
	int i = int(rnd.Rannyu(0,31)); //scelgo una delle città casualmente
	
	int indice_piu=0;
	if((i+m-1)>30){ //cioè esce a destra quindi devo usare le pbc
		indice_piu = (i+m-1)-30-1;
	}
	else{
		indice_piu = i+m-1;
	}
	
	int indice_meno=0; //non posso usare gli m indici a sinistra dell'i scelto
	if((i-m+1)<0){ //cioè esce a sinistra quindi devo usare le pbc
		indice_meno = 30 + (i-m+1) +1;
	}
	else{
		indice_meno = i-m+1;
	}	
	
	int j=0;
	if(indice_piu < indice_meno){
		do{
       			j = int(rnd.Rannyu(0,31));
     		}while((j<=indice_piu)||(j>=indice_meno));
	}
	else{	
		do{
       			j = int(rnd.Rannyu(0,31));
     		}while((j>=indice_meno)&&(j<=indice_piu));
	}
	
	vector<City> copia; 
	for (auto it = this->percorso.begin(); it != this->percorso.end(); it++) {
		copia.push_back(*it);
	}

	for(int l=0;l<m;l++){
		int pos_i = i+l;
		int pos_j = j+l;
		if(pos_i>30){
			pos_i = pos_i-30-1;
		}
		if(pos_j>30){
			pos_j = pos_j-30-1; 
		}
		swap(this->percorso[pos_i],copia[pos_j]);
		swap(this->percorso[pos_j],copia[pos_i]);
	}
	copia.clear();
	rnd.SaveSeed();
}

//Scelgo casualmente un numero m<=N-1 e inverto m città dopo la i-esima scelta casualmente.
void Individuo::mutation_inversion(Random & rnd){
	int m = int(rnd.Rannyu(0,31));
	int i = int(rnd.Rannyu(0,31)); //scelgo una delle città casualmente
	int j = i+m-1;
	if(j<=30){
		int l = int(m/2);
		for (int k=0;k<l;k++){
			swap(this->percorso[i+k],this->percorso[j-k]);
		}
	}
	else{
		j = j-30-1;
		int l = int(m/2);
		for (int k=0;k<l;k++){
			int pos_i = i+k;
			int pos_j = j-k;
			if(pos_i>30){
				pos_i = pos_i-30-1;
			}
			if(pos_j<0){
				pos_j = 30 + pos_j +1 ;
			}
			swap(this->percorso[pos_i],this->percorso[pos_j]);
		}
	}
	rnd.SaveSeed();
}

//Calcolo distanza2, tra la città i-esima e la città j-esima
double Individuo::distanza2(int i, int j) const {
	double dist2 = (this->percorso[i].getx()-this->percorso[j].getx())*(this->percorso[i].getx()-this->percorso[j].getx())+(this->percorso[i].gety()-this->percorso[j].gety())*(this->percorso[i].gety()-this->percorso[j].gety());
	return sqrt(dist2);
}

//Calcolo L1, oltre alle distanze tra le città in percorso, devo aggiungere anche la distanza tra (partenza e 0 di percorso) e tra (partenza e 30 di percorso).
double Individuo::L1() const{
	double lunghezza2=0.;
 	int n = this->percorso.size()-1;
	
  	for (int i=0; i<n; i++){
    		lunghezza2 = lunghezza2 + distanza2(i,i+1);
  	}
	
	double d_partenza_0 = sqrt((this->partenza.getx()-this->percorso[0].getx())*(this->partenza.getx()-this->percorso[0].getx())+(this->partenza.gety()-this->percorso[0].gety())*(this->partenza.gety()-this->percorso[0].gety()));

	double d_fine_partenza = sqrt((this->partenza.getx()-this->percorso[30].getx())*(this->partenza.getx()-this->percorso[30].getx())+(this->partenza.gety()-this->percorso[30].gety())*(this->partenza.gety()-this->percorso[30].gety()));
	
	lunghezza2 = lunghezza2 + d_partenza_0 + d_fine_partenza;
  
	return lunghezza2;
}

//Controllo che ogni percorso non passi mai due volte nella stessa città. Il fatto che il percorso inizi e finisca nella stessa città è implicito per come ho costruito la classe.
bool Individuo::check() {
	bool controllo = false;
	int n = this->percorso.size();
	for (int i=0; i<n-1; i++){
		for(int j=i+1; j<n; j++){
    			if((this->percorso[i].getx() == this->percorso[j].getx()) && (this->percorso[i].gety() == this->percorso[j].gety()) ){
				controllo = true;
				cout << "Il percorso passa due volte nella stessa città! Non va bene!" << endl;
			}
			else{}
		}
  	}
	return controllo;
}

//Creo il percorso totale di 32 città con in prima posizione la città di partenza (e fine).
void Individuo::set_path_totale() {
	this->path_totale.push_back(this->partenza);
	auto it1 = this->percorso.begin();
	auto it2 = this->percorso.end();
	for (auto it = it1; it !=it2; it++) {
		this->path_totale.push_back(*it);
	}
}

//Popolazione
Popolazione::Popolazione(int number, vector<Individuo>::iterator it1, vector<Individuo>::iterator it2) : N_pop(number) {
	for (auto it = it1; it !=it2; it++) {
		if((it)->check() == false){ //false per me vuol dire che non passa nella stessa città più volte
			pop.push_back(*it);
		}
		else {
			cout << "Il percorso passa due volte nella stessa città! Non va bene!" << endl;
		}
		//cout << (it)->getpercorso()[0].getx() << " " << (it)->getpercorso()[0].gety() << endl;
	}
}

Popolazione::Popolazione(int number, Random & rnd, Individuo primo) : N_pop(number) {
	if(primo.check() == false){
		pop.push_back(primo); //0
	}
	Individuo mutato(primo); //per ora è uguale a primo
	mutato.mutation_swap(rnd);
	if(mutato.check() == false){
		pop.push_back(mutato); //1
	}
	mutato.mutation_swap(rnd);
	if(mutato.check() == false){
		pop.push_back(mutato); //2
	}
	mutato.mutation_swap(rnd);
	if(mutato.check() == false){
		pop.push_back(mutato); //3
	}
	for (int i=0; i<24; i++){ //così ciclo 24 volte (da 0 a 23) e in ogni ciclo uso le 4 mutazioni. In totale creo così altre 96 mutazioni che con le altre 4 arriva a 100	
		mutato.mutation_swap(rnd);
		if(mutato.check() == false){
			pop.push_back(mutato);
		}
		mutato.mutation_shift(rnd);
		if(mutato.check() == false){
			pop.push_back(mutato);
		}
		mutato.mutation_permutation(rnd);
		if(mutato.check() == false){
			pop.push_back(mutato);
		}
		mutato.mutation_inversion(rnd);
		if(mutato.check() == false){
			pop.push_back(mutato);
		}
		
	}
	rnd.SaveSeed();
}

//Ordino gli elementi del vettore di Individui in Popolazione in base a L1
void Popolazione::sort_L1() {
	sort(this->pop.begin(), this->pop.end(), compareL1); 
}

//Selezione
int Popolazione::selezione(Random & rnd) {
	int N = this->N_pop;
	double r = rnd.Rannyu();
	double p = 3.;
	int j = int(N*pow(r,p));
	rnd.SaveSeed();
	return j;
}

//crossover: scelgo due individui della popolazione, un padre e una madre, tramite l'operatore di selezione e poi faccio il crossover
void Popolazione::crossover(Random & rnd, int k, int l) {
	Individuo madre(this->pop[k]);
	Individuo padre(this->pop[l]);

	int cut = int(rnd.Rannyu(1,30)); //scelgo una delle città casualmente, tra le 31 che ci sono nel percorso, la partenza non la considero tanto è sempre la stessa. Fino a cut escluso conservo
	vector<City> percorso_f1;
	vector<City> percorso_f2;

	//copio nei figli le prime parti fino al cut escluso
	for (int i=0; i<cut; i++) {
		percorso_f1.push_back(madre.get_percorso()[i]);
	}

	for (int i=0; i<cut; i++) {
		percorso_f2.push_back(padre.get_percorso()[i]);
	}
	
	vector<City> after_cut1;
	vector<City> after_cut2;

	//negli after cut metto tutti quelli che mancano di madre e padre
	for (int i=cut; i<31; i++) {
		after_cut1.push_back(madre.get_percorso()[i]);
		after_cut2.push_back(padre.get_percorso()[i]);
	}

	int dim = after_cut1.size(); //tanto è uguale
	
	//scorro l'altro vettore per trovare le posizioni corrispondenti
	for (int i=0; i<dim; i++) {
		for(int j=0; j<31; j++) {//devo scorrere tutto il padre per cercare
			if(after_cut1[i].getx() == padre.get_percorso()[j].getx()){
				int posiz = j;
				after_cut1[i].set_pos(posiz);
			}
			else{
			}
		}
	}
	dim = after_cut2.size();
	for (int i=0; i<dim; i++) {
		for(int j=0; j<31; j++) {//devo scorrere tutto la madre per cercare
			if(after_cut2[i].getx() == madre.get_percorso()[j].getx()){
				int posiz=j;
				after_cut2[i].set_pos(posiz);
			}
			else{
			}
		}
	}
	
	dim = after_cut2.size();
	
	vector<City> appoggio1;
	vector<City> appoggio2;

	for(int i=0; i<after_cut1.size()-1; i++){
	int min = i;
		for(int j=i+1; j<after_cut1.size(); j++){
  			if((after_cut1[j].get_pos())<(after_cut1[min].get_pos()) ){ 
    			min = j;
			}
		}

	appoggio1.push_back(after_cut1[min]);
	after_cut1.erase(after_cut1.begin()+min);
	
	i=i-1;
	}

	appoggio1.push_back(after_cut1[0]);
	auto it1 = appoggio1.begin();
	auto it2 = appoggio1.end();

	for(auto it=it1; it!=it2;it++){
		int posiz = 1;
		(*it).set_pos(posiz);
		percorso_f1.push_back(*it);
	}

	for(int i=0; i<after_cut2.size()-1; i++){
	int min = i;
		for(int j=i+1; j<after_cut2.size(); j++){
  			if((after_cut2[j].get_pos())<(after_cut2[min].get_pos()) ){
    			min = j;
			}
		}

	appoggio2.push_back(after_cut2[min]);
	after_cut2.erase(after_cut2.begin()+min);

	i=i-1;
	}
	
	appoggio2.push_back(after_cut2[0]);
	it1 = appoggio2.begin();
	it2 = appoggio2.end();
	
	for(auto it=it1; it!=it2;it++){
		int posiz = 1;
		(*it).set_pos(posiz);

		percorso_f2.push_back(*it);
	}

	dim = percorso_f1.size(); //tanto è uguale
	if (dim!=31){
		cout << "ERRORE!!" << endl;
	}
	dim = percorso_f2.size(); //tanto è uguale
	if (dim!=31){
		cout << "ERRORE!!" << endl;
	}

	//Ora sostituisco la madre e il padre con i due figli
	int getN = madre.get_N();
	City getpart = madre.get_partenza();

	it1 = percorso_f1.begin();
	it2 = percorso_f1.end();
	Individuo figlio1(getN, getpart, it1,it2);

	
	City getpart2 = padre.get_partenza();
	
	it1 = percorso_f2.begin();
	it2 = percorso_f2.end();
	Individuo figlio2(getN, getpart, it1,it2);
	
	if(figlio1.check() == false){
			this->pop.push_back(figlio1);
	}
	if(figlio2.check() == false){
			this->pop.push_back(figlio2);
	}
	
	rnd.SaveSeed();
	percorso_f1.clear();
	percorso_f2.clear();
	after_cut1.clear();
	after_cut2.clear();
	appoggio1.clear();
	appoggio2.clear();
}

//Quindi ho la prima popolazione (generazione) ordinata con tutti i metodi a posto. Ora per creare le varie nuove generazioni applico le mutazioni e il crossover con le loro rispettive probabilità
//Per implementare l'elitarismo semplicimente faccio partire le mutazioni dal secondo individuo (posizione1) ma in ogni caso il migliore (in posizione 0) può partecipare al crossover
void Popolazione::new_generazione(Random & rnd) {
	double prob_mutazione = 0.09;
	double prob_crossover = 0.6;
	double n_r = 0.;
	
	vector <Individuo> new_popol;
	for(int j=0;j<50;j++){

		int k = this->selezione(rnd);
		int l = this->selezione(rnd);
		n_r = rnd.Rannyu();
	
		//crossover
		if(n_r < prob_crossover){
			this->crossover(rnd,k,l);
			new_popol.push_back(this->pop[this->pop.size()-2]);
			new_popol.push_back(this->pop[this->pop.size()-1]);
		}
		else{
			new_popol.push_back(this->pop[l]);
			new_popol.push_back(this->pop[k]);
		}
		
		n_r = rnd.Rannyu();
		if(n_r < prob_mutazione){
			int quale = rnd.Rannyu(0,4);
   			if (quale==0) new_popol[new_popol.size()-2].mutation_swap(rnd);
      			else if (quale==1) new_popol[new_popol.size()-2].mutation_shift(rnd);
      			else if (quale==2) new_popol[new_popol.size()-2].mutation_permutation(rnd);
      			else new_popol[new_popol.size()-2].mutation_inversion(rnd);
		}

		n_r = rnd.Rannyu();
		if(n_r < prob_mutazione){
			int quale = rnd.Rannyu(0,4);
   			if (quale==0) new_popol[new_popol.size()-1].mutation_swap(rnd);
      			else if (quale==1) new_popol[new_popol.size()-1].mutation_shift(rnd);
      			else if (quale==2) new_popol[new_popol.size()-1].mutation_permutation(rnd);
      			else new_popol[new_popol.size()-1].mutation_inversion(rnd);
		}	
			
	}	
	
	this->pop.clear();
	
	for(int i=0; i<100; i++){
		this->pop.push_back(new_popol[i]);
	}

	new_popol.clear();
		
	rnd.SaveSeed();
		
}
void Popolazione::migrazione(double * new_x, double * new_y){
	City inizio(new_x[0],new_y[0]);
	vector<City> cities_nuove;
	int N=32;
	for (int i=1; i<N; i++){
		City c(new_x[i],new_y[i]);
		cities_nuove.push_back(c);
	}
	Individuo nuovo(N, inizio, cities_nuove.begin(), cities_nuove.end());
	this->pop.erase(pop.begin());
	this->pop.insert(pop.begin(),nuovo);
	cout << "fatto" << endl;
}
//Calcolo la media della metà migliore
double Popolazione::mean_L1() {
	double sum = 0.;
	for(int i=0; i<50; i++){
		sum = sum + this->pop[i].L1();
	}
	return sum/50.;
}
