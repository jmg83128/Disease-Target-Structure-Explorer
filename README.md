# 질병 표적 추천·구조 탐색기

한국어 질병명이나 KCD 코드를 건강보험심사평가원(HIRA) 질병정보 OpenAPI에서 조회해 후보를 고릅니다. API가 영문명을 반환하면 Open Targets 검색에 사용하고, 영문명이 없으면 API의 한글 표준명을 검색어로 사용합니다. 자동 기계번역은 하지 않습니다. 이후 표적 후보 재순위화, PDB 구조 비교·3D 시각화, PubMed MeSH 문헌 탐색을 제공합니다.

## Windows 실행

```bat
python -m venv .venv
.venv\Scripts\activate
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

## Google Colab 실행

Colab에서 `app.py`를 직접 실행하지 말고 Streamlit 서버로 실행합니다. `app.py`, `requirements.txt`, `kcd_terms.csv`를 Colab에 올린 뒤:

```python
!pip install -r requirements.txt
!nohup streamlit run app.py --server.port 8501 --server.address 0.0.0.0 > /content/streamlit.log 2>&1 &
```

새 셀에서 터널을 엽니다:

```python
!npx -y localtunnel --port 8501
```

## GitHub에 올리고 Streamlit으로 배포

GitHub는 코드를 저장하는 곳이고, Streamlit Community Cloud는 그 저장소의 `app.py`를 실행해 웹 주소를 만들어 주는 배포 서비스입니다. GitHub Pages만으로는 이 Python Streamlit 앱을 실행할 수 없습니다. Colab + LocalTunnel은 임시 시연용이고, 계속 열어둘 앱에는 Streamlit Community Cloud 같은 호스팅 서비스를 사용합니다.

1. `app.py`, `requirements.txt`, `README.md`, `kcd_terms.csv`를 GitHub 저장소에 올립니다.
2. [Streamlit Community Cloud](https://share.streamlit.io/)에서 GitHub 저장소와 `app.py`를 선택해 배포합니다.
3. 앱 설정의 **Secrets**에 아래 항목을 입력합니다.

```toml
HIRA_API_KEY = "공공데이터포털에서_발급받은_키"
```

앱은 이 값을 `st.secrets`에서 읽습니다. `.streamlit/secrets.toml`이나 `.env`에 실제 키를 저장하더라도 GitHub에 커밋하면 안 됩니다. `.gitignore`에 두 파일을 제외하도록 설정해 두었습니다.

배포자가 Secrets에 키를 설정하면 방문자는 자기 키를 입력할 필요가 없습니다. 대신 모든 방문자의 HIRA 조회가 배포자의 API 키와 사용량을 공유하므로, 공개 앱은 호출량 제한이나 캐시를 고려해야 합니다. 개발·시연 단계에서는 먼저 접근을 제한하고 충분히 확인한 다음 공개 URL을 공유하는 편이 안전합니다. 운영 중 트래픽 한도와 사용 조건은 [공공데이터포털 API 안내](https://www.data.go.kr/data/15119055/openapi.do)에서 확인하세요.

## 질병명 입력 방식

1. 로컬에서는 HIRA API 키를 화면에 입력하거나 `HIRA_API_KEY` 환경변수로 설정합니다. 배포 환경에서는 Streamlit Secrets에 설정합니다.
2. 한국어 질병명 또는 KCD 코드를 입력하고 `질병명·KCD 코드 찾기`를 누릅니다.
3. 반환된 코드·한글명·영문명을 확인하고 원하는 항목을 선택합니다.
4. API의 영문명(없으면 한글명)을 Open Targets 검색어로 넣고, Open Targets 질병 후보 중 정확한 표준 ID를 선택합니다.

HIRA API에 연결할 수 없을 때를 위해 CSV 용어집도 보조로 지원합니다. CSV는 `kcd_code`, `korean_name`, `english_name` 열이 필요하고, 선택적으로 `aliases` 열을 둘 수 있습니다. `aliases`에는 한국어 검색어를 세미콜론(`;`)으로 구분해 적습니다. `kcd_terms.csv`를 `app.py`와 같은 폴더에 두면 자동으로 읽고, 화면에서 CSV를 추가 업로드할 수도 있습니다. 현재 L20.8 행은 대화에서 제시된 예시이며 공식 KCD 원본 전체와 대조한 자료가 아닙니다. L20.8은 기타 아토피성 피부염의 세부 코드이므로 모든 아토피 질환을 뜻하지 않습니다. 앱은 후보를 보여줄 뿐 해당 코드가 사용자의 진단에 맞는다고 판정하지 않습니다.

영어명을 직접 입력해 검색할 수도 있습니다. 로컬 환경변수로 설정하려면 Windows 명령 프롬프트에서 `set HIRA_API_KEY=발급받은키`를 실행한 뒤 앱을 실행하세요. API 사용승인 상태와 호출 한도는 공공데이터포털에서 확인해야 합니다.

## 기타 기능

- Open Targets 상위 25개 후보를 질병 연관 점수, genetics 점수, known drug 점수, tractability 신호의 사용자 가중 합산으로 재정렬합니다.
- 표적 후보의 UniProt accession으로 RCSB PDB 실험 구조를 찾고 실험 방법, 해상도, mutation 주석, 기준 서열 대비 차이 후보·커버리지 및 비단백질 성분을 보여줍니다.
- 선택 구조를 회전 가능한 3D 뷰어로 표시합니다.
- 질병명으로 PubMed MeSH 필드를 검색하고 논문 링크·PMID·DOI를 표시합니다.

## 해석상의 제한

- 표적 재순위 점수는 설명 가능한 가중 합산 기준선이며 학습된 머신러닝 모델이 아닙니다. Open Targets 상위 25개로 입력 후보가 제한됩니다.
- Open Targets 연관 점수는 확률이나 인과성 점수가 아닙니다. 근거 자료의 공백도 점수에 반영될 수 있습니다.
- 기준 서열 차이 후보는 isoform, 종 차이, construct 설계일 수 있어 실제 mutation으로 단정하지 않습니다.
- 해상도는 X-ray와 Cryo-EM 등 방법 간 단순 비교하지 않습니다.
- 비단백질 성분이 모두 약물 리간드는 아닙니다. 자동 검색만으로 도킹 pocket 적합성을 확정하지 않습니다.
- MeSH 문헌 검색은 질병 ID를 MeSH ID로 자동 변환하는 기능이 아닙니다.

## 데이터 출처

- 질병명·코드 검색: [건강보험심사평가원 질병정보서비스 OpenAPI](https://www.data.go.kr/data/15119055/openapi.do) (KCD 기반 질병정보; 공공데이터포털 서비스 키 필요)
- Open Targets GraphQL: [공식 문서](https://platform-docs.opentargets.org/data-access/graphql-api)
- RCSB PDB Search API: [공식 문서](https://search.rcsb.org/)
- RCSB PDB Data API: [공식 문서](https://data.rcsb.org/)
- UniProt REST API: [공식 문서](https://www.uniprot.org/help/api)
- NCBI Entrez E-utilities: [공식 문서](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
