#ifndef AnalyzeIfValid_h
#define AnalyzeIfValid_h

template <class A, class T1>
void analyzeIfValid(A & analyser, const edm::Handle<T1> & h1, TTree * tree) {
  if (h1.isValid())
    analyser.analyze(*h1, tree);
}

template <class A, class T1, class T2>
void analyzeIfValid(A & analyser, const edm::Handle<T1> & h1, const edm::Handle<T2> & h2, TTree * tree) {
  if (h1.isValid() and h2.isValid())
    analyser.analyze(*h1, *h2, tree);
}

template <class A, class T1, class T2, class T3>
void analyzeIfValid(A & analyser, const edm::Handle<T1> & h1, const edm::Handle<T2> & h2, const edm::Handle<T3> & h3, TTree * tree) {
  if (h1.isValid() and h2.isValid() and h3.isValid())
    analyser.analyze(*h1, *h2, *h3, tree);
}

template <class A, class T1, class T2, class T3, class T4>
void analyzeIfValid(A & analyser, const edm::Handle<T1> & h1, const edm::Handle<T2> & h2, const edm::Handle<T3> & h3, const edm::Handle<T4> & h4, TTree * tree) {
  if (h1.isValid() and h2.isValid() and h3.isValid() and h4.isValid())
    analyser.analyze(*h1, *h2, *h3, *h4, tree);
}

template <class A, class T1, class T2, class T3, class T4, class T5>
void analyzeIfValid(A & analyser, const edm::Handle<T1> & h1, const edm::Handle<T2> & h2, const edm::Handle<T3> & h3, const edm::Handle<T4> & h4, const edm::Handle<T5> & h5, TTree * tree) {
  if (h1.isValid() and h2.isValid() and h3.isValid() and h4.isValid() and h5.isValid())
    analyser.analyze(*h1, *h2, *h3, *h4, *h5, tree);
}

template <class A, class T1, class T2, class T3, class T4, class T5, class T6>
void analyzeIfValid(A & analyser, const edm::Handle<T1> & h1, const edm::Handle<T2> & h2, const edm::Handle<T3> & h3, const edm::Handle<T4> & h4, const edm::Handle<T5> & h5, const edm::Handle<T6> & h6, TTree * tree) {
  if (h1.isValid() and h2.isValid() and h3.isValid() and h4.isValid() and h5.isValid() and h6.isValid())
    analyser.analyze(*h1, *h2, *h3, *h4, *h5, *h6, tree);
}

template <class A, class T1, class T2, class T3, class T4, class T5, class T6, class T7>
void analyzeIfValid(A & analyser, const edm::Handle<T1> & h1, const edm::Handle<T2> & h2, const edm::Handle<T3> & h3, const edm::Handle<T4> & h4, const edm::Handle<T5> & h5, const edm::Handle<T6> & h6, const edm::Handle<T7> & h7, TTree * tree) {
  if (h1.isValid() and h2.isValid() and h3.isValid() and h4.isValid() and h5.isValid() and h6.isValid() and h7.isValid())
    analyser.analyze(*h1, *h2, *h3, *h4, *h5, *h6, *h7, tree);
}

template <class A, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8>
void analyzeIfValid(A & analyser, const edm::Handle<T1> & h1, const edm::Handle<T2> & h2, const edm::Handle<T3> & h3, const edm::Handle<T4> & h4, const edm::Handle<T5> & h5, const edm::Handle<T6> & h6, const edm::Handle<T7> & h7, const edm::Handle<T8> & h8, TTree * tree) {
  if (h1.isValid() and h2.isValid() and h3.isValid() and h4.isValid() and h5.isValid() and h6.isValid() and h7.isValid() and h8.isValid())
    analyser.analyze(*h1, *h2, *h3, *h4, *h5, *h6, *h7, *h8, tree);
}

template <class A, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
void analyzeIfValid(A & analyser, const edm::Handle<T1> & h1, const edm::Handle<T2> & h2, const edm::Handle<T3> & h3, const edm::Handle<T4> & h4, const edm::Handle<T5> & h5, const edm::Handle<T6> & h6, const edm::Handle<T7> & h7, const edm::Handle<T8> & h8, const edm::Handle<T9> & h9, TTree * tree) {
  if (h1.isValid() and h2.isValid() and h3.isValid() and h4.isValid() and h5.isValid() and h6.isValid() and h7.isValid() and h8.isValid() and h9.isValid())
    analyser.analyze(*h1, *h2, *h3, *h4, *h5, *h6, *h7, *h8, *h9, tree);
}

template <class A, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9, class T10>
void analyzeIfValid(A & analyser, const edm::Handle<T1> & h1, const edm::Handle<T2> & h2, const edm::Handle<T3> & h3, const edm::Handle<T4> & h4, const edm::Handle<T5> & h5, const edm::Handle<T6> & h6, const edm::Handle<T7> & h7, const edm::Handle<T8> & h8, const edm::Handle<T9> & h9, const edm::Handle<T10> & h10, TTree * tree) {
  if (h1.isValid() and h2.isValid() and h3.isValid() and h4.isValid() and h5.isValid() and h6.isValid() and h7.isValid() and h8.isValid() and h9.isValid() and h10.isValid())
    analyser.analyze(*h1, *h2, *h3, *h4, *h5, *h6, *h7, *h8, *h9, *h10, tree);
}

#endif // AnalyzeIfValid_h
