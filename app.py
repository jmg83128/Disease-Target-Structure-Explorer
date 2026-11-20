import csv
import os
import re
import xml.etree.ElementTree as ET
from io import StringIO
from pathlib import Path
from urllib.parse import unquote

import py3Dmol
import requests
import streamlit as st
import streamlit.components.v1 as components
from Bio.Align import PairwiseAligner


st.set_page_config(page_title="질병 표적 추천·구조 탐색기", page_icon="🧬", layout="wide")
st.title("🧬 질병 표적 추천·구조 탐색기")
st.caption("질병 근거로 표적 후보를 재순위화하고, PDB 구조를 비교·시각화하는 연구용 프로토타입")

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
UNIPROT_REST = "https://rest.uniprot.org/uniprotkb"
HIRA_DISEASE_API = "https://apis.data.go.kr/B551182/diseaseInfoService1/getDissNameCodeList1"
RCSB_SEARCH_ATTRIBUTE = "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers"
BUILTIN_KCD = [
    {
        "kcd_code": "L20.8",
        "korean_name": "기타 아토피성 피부염",
        "english_name": "Other atopic dermatitis",
        "aliases": "아토피성 피부염",
        "source": "사용자 제공 검색 예시 (공식 KCD 원본 대조 필요)",
    }
]


def contains_korean(text: str) -> bool:
    return bool(re.search(r"[가-힣ㄱ-ㅎㅏ-ㅣ]", text))


def configured_hira_key() -> str:
    """Read the deployer's key from Streamlit secrets or a local env var."""
    try:
        secret_key = st.secrets.get("HIRA_API_KEY", "")
    except Exception:
        secret_key = ""
    return str(secret_key or os.getenv("HIRA_API_KEY", "")).strip()


def is_kcd_code(text: str) -> bool:
    return bool(re.fullmatch(r"[A-Za-z][0-9]{2}(?:\.[0-9A-Za-z]{1,2})?", text.strip()))


def parse_kcd_text(raw: bytes, source: str) -> list[dict]:
    try:
        content = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        content = raw.decode("cp949")
    reader = csv.DictReader(StringIO(content))
    if not reader.fieldnames:
        raise ValueError("CSV 헤더를 읽을 수 없습니다.")

    header_lookup = {header.strip().lower(): header for header in reader.fieldnames if header}

    def column(*options):
        return next((header_lookup[name.lower()] for name in options if name.lower() in header_lookup), None)

    code_col = column("kcd_code", "kcd", "code", "코드")
    korean_col = column("korean_name", "korean", "한글명", "국문명")
    english_col = column("english_name", "english", "영문명")
    alias_col = column("aliases", "korean_aliases", "synonyms", "동의어", "별칭")
    if not code_col or not korean_col or not english_col:
        raise ValueError("CSV에 kcd_code, korean_name, english_name 열이 필요합니다. aliases 열은 선택입니다.")

    records = []
    for row in reader:
        code = (row.get(code_col) or "").strip().upper()
        korean = (row.get(korean_col) or "").strip()
        english = (row.get(english_col) or "").strip()
        if not (code and korean and english):
            continue
        records.append(
            {
                "kcd_code": code,
                "korean_name": korean,
                "english_name": english,
                "aliases": (row.get(alias_col) or "").strip() if alias_col else "",
                "source": source,
            }
        )
    return records


def parse_kcd_csv(uploaded_file) -> list[dict]:
    return parse_kcd_text(uploaded_file.getvalue(), "업로드한 KCD CSV")


def find_kcd_matches(query: str, records: list[dict]) -> list[dict]:
    normalized_query = re.sub(r"\s+", " ", query.strip()).casefold()
    if not normalized_query:
        return []
    exact, partial = [], []
    for record in records:
        names = [record.get("korean_name", ""), record.get("english_name", "")]
        aliases = re.split(r"[;|,，、]", record.get("aliases", ""))
        names.extend(aliases)
        code = record.get("kcd_code", "").strip().casefold()
        normalized_names = [re.sub(r"\s+", " ", name.strip()).casefold() for name in names if name.strip()]
        if normalized_query == code or normalized_query in normalized_names:
            exact.append(record)
        elif any(normalized_query in name or name in normalized_query for name in normalized_names):
            partial.append(record)
    return exact + partial


@st.cache_data(ttl=86400, show_spinner=False)
def search_hira_diseases(query: str, service_key: str, limit: int = 30) -> list[dict]:
    """Search HIRA's KCD-based disease information OpenAPI."""
    params = {
        # data.go.kr often provides an already URL-encoded key. Decode once
        # before requests encodes query parameters, avoiding %252F/%252B.
        "serviceKey": unquote(service_key.strip()),
        "sickType": "2",  # four-character/most detailed level
        "medTp": "1",  # western medicine
        "diseaseType": "SICK_NM" if contains_korean(query) else "SICK_CD",
        "searchText": query.strip(),
        "numOfRows": limit,
        "pageNo": 1,
    }
    response = requests.get(HIRA_DISEASE_API, params=params, timeout=30)
    if not response.ok:
        # requests.HTTPError includes the complete request URL (and service key).
        # Never surface that URL in the Streamlit error message.
        raise RuntimeError(
            f"HIRA API가 HTTP {response.status_code}를 반환했습니다. "
            "서비스 키 인코딩, API 활용 신청/승인 상태, 해당 API의 접근 권한을 확인해 주세요. "
            "보안을 위해 요청 URL과 키는 오류 메시지에 표시하지 않았습니다."
        )
    content_type = response.headers.get("Content-Type", "").lower()
    if "json" in content_type or response.text.lstrip().startswith("{"):
        payload = response.json()
        api_response = payload.get("response", {})
        header = api_response.get("header", {})
        items = api_response.get("body", {}).get("items", {}).get("item", [])
    else:
        root = ET.fromstring(response.content)

        def local_name(element):
            return element.tag.rsplit("}", 1)[-1]

        header_node = next((node for node in root.iter() if local_name(node) == "header"), None)
        header = {local_name(child): child.text or "" for child in header_node} if header_node is not None else {}
        item_nodes = [node for node in root.iter() if local_name(node) == "item"]
        items = [
            {local_name(child): child.text or "" for child in node}
            for node in item_nodes
        ]
    if str(header.get("resultCode", "00")) not in {"00", "0000"}:
        raise RuntimeError(header.get("resultMsg") or "HIRA API 요청에 실패했습니다.")
    if isinstance(items, dict):
        items = [items]
    results = []
    for item in items or []:
        # HIRA versions may return a list of ICD/KCD descriptors. Keep the
        # original fields for debugging while normalizing common field names.
        code = str(item.get("sickCd") or item.get("SICK_CD") or item.get("diseaseCode") or "").strip()
        korean = str(item.get("sickNm") or item.get("SICK_NM") or item.get("diseaseName") or "").strip()
        english = str(
            item.get("sickEngNm") or item.get("SICK_ENG_NM") or item.get("engSickNm")
            or item.get("ENG_SICK_NM") or item.get("englishName") or ""
        ).strip()
        if code or korean or english:
            results.append({
                "kcd_code": code,
                "korean_name": korean,
                "english_name": english,
                "source": "건강보험심사평가원 질병정보 OpenAPI (KCD 기반)",
                "raw": item,
            })
    return results


def graphql_post(url: str, query: str, variables: dict) -> dict:
    response = requests.post(url, json={"query": query, "variables": variables}, timeout=60)
    response.raise_for_status()
    payload = response.json()
    if payload.get("errors"):
        raise RuntimeError(payload["errors"][0].get("message", "GraphQL query failed"))
    return payload.get("data", {})


@st.cache_data(ttl=86400, show_spinner=False)
def search_diseases(query_text: str) -> list[dict]:
    query = """
    query SearchDisease($queryString: String!) {
      search(queryString: $queryString, entityNames: ["disease"]) {
        hits { id name entity description }
      }
    }
    """
    # HIRA code labels may include modifiers that are absent from EFO labels,
    # e.g. "Eczema herpeticum" vs the broader "eczema" concept. Search the
    # exact API term first, then clearly marked broader lexical forms.
    search_terms = [query_text.strip()]
    leading_qualifier_removed = re.sub(
        r"^(?:other|unspecified|unclassified|not elsewhere classified|not otherwise specified)\s+",
        "",
        query_text.strip(),
        flags=re.IGNORECASE,
    ).strip()
    modifier_removed = re.sub(
        r"\b(?:herpes\s+simplex|herpes|herpeticum|herpetic|simplex)\b",
        " ",
        query_text.strip(),
        flags=re.IGNORECASE,
    )
    modifier_removed = re.sub(r"\s+", " ", modifier_removed).strip(" ,;-")
    for broader in (leading_qualifier_removed, modifier_removed):
        if broader and broader.casefold() not in {term.casefold() for term in search_terms}:
            search_terms.append(broader)

    combined = []
    seen_ids = set()
    for term in search_terms:
        data = graphql_post(OT_GRAPHQL, query, {"queryString": term})
        for hit in data.get("search", {}).get("hits", []):
            hit_id = hit.get("id")
            if not hit_id or hit_id in seen_ids:
                continue
            seen_ids.add(hit_id)
            hit["matched_query"] = term
            hit["is_broader_match"] = term.casefold() != query_text.strip().casefold()
            combined.append(hit)
    return combined


def get_targets(disease_id: str, page_size: int = 25) -> dict:
    query = """
    query DiseaseTargets($diseaseId: String!, $size: Int!) {
      disease(efoId: $diseaseId) {
        id
        name
        associatedTargets(page: {index: 0, size: $size}, orderByScore: "score") {
          count
          rows {
            target { id approvedSymbol approvedName proteinIds { id source } tractability { label modality value } }
            score
            datatypeScores { id score }
          }
        }
      }
    }
    """
    data = graphql_post(OT_GRAPHQL, query, {"diseaseId": disease_id, "size": page_size})
    return data.get("disease") or {
        "id": disease_id,
        "name": "",
        "associatedTargets": {"count": 0, "rows": []},
    }


def component_score(row: dict, category: str) -> float:
    aliases = {
        "genetics": ("genetic", "gwas", "gene_burden", "ot_genetics", "variant"),
        "known_drug": ("known_drug", "chembl"),
    }[category]
    values = [
        float(item.get("score") or 0)
        for item in row.get("datatypeScores", [])
        if any(alias in str(item.get("id", "")).lower() for alias in aliases)
    ]
    return max(values, default=0.0)


def has_tractability(row: dict) -> bool:
    return any(
        item.get("value") is True or str(item.get("value", "")).lower() == "true"
        for item in row.get("target", {}).get("tractability", [])
    )


def rank_targets(rows: list[dict], weights: dict[str, int]) -> list[dict]:
    weight_sum = sum(weights.values()) or 1
    ranked = []
    for row in rows:
        genetic = component_score(row, "genetics")
        known_drug = component_score(row, "known_drug")
        tractability = 1.0 if has_tractability(row) else 0.0
        association = float(row.get("score") or 0)
        priority = (
            weights["association"] * association
            + weights["genetics"] * genetic
            + weights["known_drug"] * known_drug
            + weights["tractability"] * tractability
        ) / weight_sum
        target = row.get("target", {})
        protein_ids = target.get("proteinIds") or []
        uniprot = next((p.get("id") for p in protein_ids if p.get("source") == "uniprot_swissprot"), None)
        if not uniprot:
            uniprot = next((p.get("id") for p in protein_ids if str(p.get("source", "")).startswith("uniprot")), None)
        ranked.append(
            {
                "priority": priority,
                "association": association,
                "genetics": genetic,
                "known_drug": known_drug,
                "tractability": tractability,
                "target": target,
                "uniprot": uniprot,
            }
        )
    return sorted(ranked, key=lambda item: item["priority"], reverse=True)


def search_pubmed_mesh(disease_name: str, limit: int = 8) -> list[dict]:
    response = requests.get(
        f"{NCBI_EUTILS}/esearch.fcgi",
        params={
            "db": "pubmed",
            "term": f'"{disease_name}"[MeSH Terms]',
            "retmode": "json",
            "retmax": limit,
            "sort": "relevance",
        },
        timeout=30,
    )
    response.raise_for_status()
    ids = response.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    summary = requests.get(
        f"{NCBI_EUTILS}/esummary.fcgi",
        params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"},
        timeout=30,
    )
    summary.raise_for_status()
    result = summary.json().get("result", {})
    articles = []
    for pmid in ids:
        item = result.get(pmid, {})
        if not item:
            continue
        doi = next(
            (a.get("value") for a in item.get("articleids", []) if a.get("idtype") == "doi"),
            "",
        )
        articles.append(
            {
                "title": item.get("title", "제목 없음").rstrip("."),
                "authors": ", ".join(a.get("name", "") for a in item.get("authors", [])[:3]),
                "journal": item.get("fulljournalname", ""),
                "year": item.get("pubdate", "")[:4],
                "pmid": pmid,
                "doi": doi,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            }
        )
    return articles


@st.cache_data(ttl=86400, show_spinner=False)
def get_uniprot_sequence(accession: str) -> str:
    response = requests.get(f"{UNIPROT_REST}/{accession}.fasta", timeout=30)
    response.raise_for_status()
    return "".join(line.strip() for line in response.text.splitlines() if not line.startswith(">"))


def sequence_comparison(reference: str, construct: str) -> tuple[float | None, list[str]]:
    reference = re.sub(r"[^A-Z]", "", reference.upper())
    construct = re.sub(r"[^A-Z]", "", construct.upper())
    if not reference or not construct:
        return None, []
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference, construct)[0]
    coordinates = alignment.coordinates
    aligned_reference_positions = set()
    differences = []
    for block in range(coordinates.shape[1] - 1):
        r0, r1 = int(coordinates[0, block]), int(coordinates[0, block + 1])
        c0, c1 = int(coordinates[1, block]), int(coordinates[1, block + 1])
        if r1 <= r0 or c1 <= c0:
            continue
        aligned_length = min(r1 - r0, c1 - c0)
        for offset in range(aligned_length):
            rpos, cpos = r0 + offset, c0 + offset
            aligned_reference_positions.add(rpos)
            if reference[rpos] != construct[cpos] and reference[rpos] != "X" and construct[cpos] != "X":
                differences.append(f"{reference[rpos]}{rpos + 1}{construct[cpos]}")
    coverage = len(aligned_reference_positions) / len(reference)
    return coverage, differences


def search_pdb_by_uniprot(accession: str, limit: int = 40) -> list[dict]:
    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": f"{RCSB_SEARCH_ATTRIBUTE}.database_accession",
                        "operator": "exact_match",
                        "value": accession,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": f"{RCSB_SEARCH_ATTRIBUTE}.database_name",
                        "operator": "exact_match",
                        "value": "UniProt",
                    },
                },
            ],
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
            "results_content_type": ["experimental"],
        },
    }
    response = requests.post(RCSB_SEARCH, json=payload, timeout=45)
    if response.status_code == 204:
        return []
    response.raise_for_status()
    hits = response.json().get("result_set", [])
    entities = []
    for hit in hits:
        identifier = hit.get("identifier", "")
        match = re.match(r"^([A-Za-z0-9]{4})_(\d+)$", identifier)
        if match:
            entities.append({"pdb_id": match.group(1).upper(), "entity_id": match.group(2)})
    return entities


def get_pdb_metadata(pdb_entities: list[dict], accession: str) -> list[dict]:
    pdb_ids = list(dict.fromkeys(item["pdb_id"] for item in pdb_entities))
    if not pdb_ids:
        return []
    query = """
    query PdbEntries($entryIds: [String!]!) {
      entries(entry_ids: $entryIds) {
        rcsb_id
        struct { title }
        exptl { method }
        rcsb_entry_info { resolution_combined }
        polymer_entities {
          rcsb_id
          entity_poly { pdbx_seq_one_letter_code_can rcsb_sample_sequence_length }
          rcsb_polymer_entity { pdbx_description pdbx_mutation pdbx_fragment }
          rcsb_polymer_entity_container_identifiers {
            entity_id
            reference_sequence_identifiers { database_name database_accession }
          }
        }
        nonpolymer_entities {
          nonpolymer_comp { chem_comp { id name } }
        }
      }
    }
    """
    data = graphql_post(RCSB_GRAPHQL, query, {"entryIds": pdb_ids})
    entries = data.get("entries", [])
    requested_entities = {(i["pdb_id"], str(i["entity_id"])) for i in pdb_entities}
    found = []
    for entry in entries:
        if not entry:
            continue
        pdb_id = entry.get("rcsb_id", "")
        for entity in entry.get("polymer_entities", []) or []:
            ids = entity.get("rcsb_polymer_entity_container_identifiers", {}) or {}
            entity_id = str(ids.get("entity_id", ""))
            if (pdb_id, entity_id) not in requested_entities:
                continue
            references = ids.get("reference_sequence_identifiers", []) or []
            matching = any(
                ref.get("database_name", "").lower() == "uniprot"
                and ref.get("database_accession", "") == accession
                for ref in references
            )
            if not matching:
                continue
            entity_poly = entity.get("entity_poly", {}) or {}
            entity_metadata = entity.get("rcsb_polymer_entity", {}) or {}
            sequence = re.sub(r"[^A-Z]", "", (entity_poly.get("pdbx_seq_one_letter_code_can") or "").upper())
            found.append(
                {
                    "pdb_id": pdb_id,
                    "entity_id": entity_id,
                    "title": (entry.get("struct") or {}).get("title", ""),
                    "methods": [x.get("method", "") for x in (entry.get("exptl") or []) if x.get("method")],
                    "resolutions": (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or [],
                    "description": entity_metadata.get("pdbx_description", ""),
                    "mutation_annotation": entity_metadata.get("pdbx_mutation") or "Not reported",
                    "fragment": entity_metadata.get("pdbx_fragment", ""),
                    "sequence": sequence,
                    "sequence_length": entity_poly.get("rcsb_sample_sequence_length") or len(sequence),
                    "components": [
                        comp.get("chem_comp", {}).get("id", "")
                        for nonpolymer in (entry.get("nonpolymer_entities") or [])
                        for comp in [nonpolymer.get("nonpolymer_comp") or {}]
                        if comp.get("chem_comp", {}).get("id")
                    ],
                }
            )
    return found


def attach_sequence_metrics(structures: list[dict], uniprot_sequence: str) -> list[dict]:
    for structure in structures:
        coverage, mismatches = sequence_comparison(uniprot_sequence, structure.get("sequence", ""))
        structure["reference_coverage"] = coverage
        structure["sequence_differences"] = mismatches
        resolution_values = structure.get("resolutions", [])
        structure["best_resolution"] = min(resolution_values) if resolution_values else None
        methods = " ".join(structure.get("methods", [])).lower()
        if "x-ray" in methods or "xray" in methods:
            method_type = "X-ray"
        elif "electron microscopy" in methods or "electron cryo" in methods or "cryo" in methods:
            method_type = "Cryo-EM"
        elif "nmr" in methods:
            method_type = "NMR"
        else:
            method_type = "; ".join(structure.get("methods", [])) or "Unknown"
        structure["method_type"] = method_type
        nonwater = [
            code for code in structure.get("components", [])
            if code not in {"HOH", "DOD", "GOL", "EDO", "PEG", "SO4", "PO4", "CL", "NA", "MG", "ZN", "CA", "K"}
        ]
        structure["nonwater_components"] = list(dict.fromkeys(nonwater))
    return structures


def render_structure(pdb_id: str):
    response = requests.get(f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=60)
    response.raise_for_status()
    viewer = py3Dmol.view(width=1000, height=580)
    viewer.addModel(response.text, "cif")
    viewer.setStyle({"cartoon": {"color": "spectrum"}})
    viewer.addStyle({"hetflag": True}, {"stick": {"radius": 0.18}})
    viewer.setBackgroundColor("0xf7f9fc")
    viewer.zoomTo()
    components.html(viewer._make_html(), height=600, scrolling=False)


def clear_after_new_search():
    for key in ("disease_hits", "target_result", "ranked_targets", "selected_disease_id", "pubmed_articles", "pdb_structures", "pdb_target_id"):
        st.session_state.pop(key, None)


if "english_query" not in st.session_state:
    st.session_state.english_query = ""
deployed_hira_key = configured_hira_key()
if "hira_service_key" not in st.session_state:
    st.session_state.hira_service_key = deployed_hira_key

st.markdown("### 질병명 표준화")
st.write("한국어명이나 KCD 코드를 HIRA 질병정보 API에서 조회합니다. API가 돌려준 영문명이 있으면 그 이름을 사용하고, 영문명이 없으면 한글 표준명을 검색어로 사용합니다. 기계번역은 하지 않습니다.")
if deployed_hira_key:
    st.session_state.hira_service_key = deployed_hira_key
    st.caption("HIRA API 키가 서버 비밀 설정에서 안전하게 연결되어 있습니다.")
else:
    st.text_input(
        "HIRA 질병정보 API 서비스 키",
        type="password",
        key="hira_service_key",
        help="로컬에서 시험할 때만 입력하세요. 공개 배포에서는 앱의 Secrets 설정을 사용하세요.",
    )
uploaded_kcd = st.file_uploader("KCD 용어 CSV 업로드 (선택)", type=["csv"], help="필수 열: kcd_code, korean_name, english_name. aliases 열에는 동의어를 세미콜론으로 구분해 넣을 수 있습니다.")
kcd_records = BUILTIN_KCD.copy()
bundled_kcd_path = Path(__file__).with_name("kcd_terms.csv")
if bundled_kcd_path.exists():
    try:
        bundled_records = parse_kcd_text(bundled_kcd_path.read_bytes(), "프로젝트 CSV 예시 (공식 KCD 원본 대조 필요)")
        bundled_codes = {row["kcd_code"] for row in bundled_records}
        kcd_records = [row for row in kcd_records if row["kcd_code"] not in bundled_codes] + bundled_records
    except Exception as exc:
        st.warning(f"프로젝트 KCD CSV를 읽지 못해 기본 예시만 사용합니다: {exc}")
if uploaded_kcd is not None:
    try:
        uploaded_records = parse_kcd_csv(uploaded_kcd)
        # Uploaded official/local data takes precedence over the example row with the same code.
        uploaded_codes = {row["kcd_code"] for row in uploaded_records}
        kcd_records = [row for row in kcd_records if row["kcd_code"] not in uploaded_codes] + uploaded_records
        st.caption(f"업로드한 KCD 항목 {len(uploaded_records):,}개를 검색에 사용합니다.")
    except Exception as exc:
        st.error(f"KCD CSV를 읽지 못했습니다: {exc}")

with st.form("disease_form"):
    disease_input = st.text_input("한국어 질병명, KCD 코드 또는 영어명", placeholder="예: 아토피 / L20.8 / Other atopic dermatitis")
    submit_lookup = st.form_submit_button("질병명·KCD 코드 찾기", type="primary")

if submit_lookup:
    clear_after_new_search()
    st.session_state.pop("kcd_matches", None)
    st.session_state.pop("kcd_source_label", None)
    st.session_state.english_query = ""
    query = disease_input.strip()
    if not query:
        st.warning("질병명이나 KCD 코드를 입력해 주세요.")
    elif contains_korean(query) or is_kcd_code(query):
        matches = []
        if st.session_state.hira_service_key.strip():
            try:
                with st.spinner("HIRA 질병정보 API에서 질병명과 코드를 찾는 중입니다..."):
                    matches = search_hira_diseases(query, st.session_state.hira_service_key)
                if matches:
                    st.session_state.kcd_source_label = "건강보험심사평가원 질병정보 OpenAPI (KCD 기반)"
            except Exception as exc:
                st.warning(f"HIRA API 조회 실패: {exc} · 로컬 CSV 용어집으로 검색합니다.")
        if not matches:
            matches = find_kcd_matches(query, kcd_records)
            if matches:
                st.session_state.kcd_source_label = "로컬 KCD CSV 용어집"
        if matches:
            st.session_state.kcd_matches = matches
        else:
            st.warning("질병 항목을 찾지 못했어요. HIRA API 서비스 키와 활용 신청 상태를 확인하거나 KCD CSV를 업로드해 주세요. 자동 번역은 하지 않습니다.")
    else:
        # English input can be searched directly, without pretending it came from a KCD mapping.
        st.session_state.english_query = query
        st.session_state.kcd_source_label = "직접 입력한 영어 검색어"

kcd_matches = st.session_state.get("kcd_matches", [])
if kcd_matches:
    st.subheader("질병명·코드 후보 확인")
    st.caption("입력 표현이 넓은 질병명이나 별칭이면 더 구체적인 코드가 후보로 보일 수 있습니다. 코드를 확인해 원하는 항목을 선택하세요.")
    kcd_options = {
        f"{row.get('kcd_code') or '코드 없음'} | {row.get('korean_name') or '한글명 없음'} | {row.get('english_name') or '영문명 없음'}": row
        for row in kcd_matches
    }
    selected_kcd_label = st.selectbox("KCD 항목", list(kcd_options))
    selected_kcd = kcd_options[selected_kcd_label]
    st.caption(f"자료: {selected_kcd.get('source', st.session_state.get('kcd_source_label', 'HIRA 질병정보 API'))}")
    if not selected_kcd.get("english_name"):
        st.info("이 API 응답에 영문명이 없어요. 한글 표준명으로 Open Targets에서 검색하거나, 영문명을 확인해 직접 입력할 수 있습니다.")
    if st.button("선택한 질병명으로 검색어 설정"):
        clear_after_new_search()
        st.session_state.english_query = selected_kcd.get("english_name") or selected_kcd.get("korean_name", "")
        st.session_state.kcd_source_label = f"{selected_kcd.get('kcd_code', '')} · {selected_kcd.get('korean_name', '')} (HIRA 질병정보 API)"

st.text_input("Open Targets 검색어 (필요하면 수정)", key="english_query")
if st.session_state.get("kcd_source_label"):
    st.caption(f"검색어 출처: {st.session_state.kcd_source_label}")
if st.button("Open Targets 질병 후보 찾기", disabled=not st.session_state.english_query.strip()):
    clear_after_new_search()
    try:
        with st.spinner("질병 후보를 찾는 중입니다..."):
            st.session_state.disease_hits = search_diseases(st.session_state.english_query.strip())
    except Exception as exc:
        st.error(f"Open Targets 검색 오류: {exc}")

hits = st.session_state.get("disease_hits", [])
selected = None
if hits:
    st.subheader("질병 개념 확인")
    st.caption("후보에 표시된 검색어도 확인하세요. 수식어를 뺀 넓은 검색어로 찾은 후보는 원래 KCD 질병과 범위가 다를 수 있으니, 맞는 질병인지 직접 확인해야 합니다.")
    hit_options = {
        f"{hit.get('name', '이름 없음')} | {hit.get('id', '')} | 검색어: {hit.get('matched_query', '')}": hit
        for hit in hits if hit.get("id")
    }
    selected_label = st.selectbox("Open Targets 질병 후보", list(hit_options))
    selected = hit_options[selected_label]
    if selected.get("is_broader_match"):
        st.warning(f"'{selected.get('matched_query')}'로 넓혀 검색해 찾은 후보입니다. 입력한 KCD 질환보다 범위가 넓거나 다를 수 있으니 확인하세요.")
    if selected.get("description"):
        st.write(selected["description"])
    if st.button("이 질병의 표적 후보 계산"):
        try:
            with st.spinner("질병-표적 근거를 불러옵니다..."):
                st.session_state.target_result = get_targets(selected["id"])
            st.session_state.selected_disease_id = selected["id"]
            st.session_state.pop("pubmed_articles", None)
            st.session_state.pop("pdb_structures", None)
        except Exception as exc:
            st.error(f"표적 조회 오류: {exc}")
elif st.session_state.get("disease_hits") == []:
    st.info("Open Targets에서 이 표현과 완화 검색어에 해당하는 질병 개념을 찾지 못했습니다. 검색어를 수정해 다시 시도해 주세요.")

target_result = st.session_state.get("target_result")
if target_result and st.session_state.get("selected_disease_id") == (selected or {}).get("id"):
    disease_id = st.session_state.selected_disease_id
    rows = target_result.get("associatedTargets", {}).get("rows", [])
    st.header(f"표적 후보 추천: {target_result.get('name') or selected.get('name')}")
    st.caption("Open Targets 상위 후보를 가져온 뒤, 아래 근거 가중치로 다시 정렬합니다. 점수는 연구용 우선순위 지표이지 치료 효과의 확률이 아닙니다.")

    with st.expander("추천 기준 조정", expanded=True):
        st.write("가중치 합은 자동으로 정규화됩니다. 근거가 제공되지 않은 항목은 0으로 계산되므로 표의 구성 점수도 함께 확인하세요.")
        w1, w2, w3, w4 = st.columns(4)
        weights = {
            "association": w1.slider("질병 연관 근거", 0, 100, 50),
            "genetics": w2.slider("유전학 근거", 0, 100, 25),
            "known_drug": w3.slider("기존 약물 근거", 0, 100, 15),
            "tractability": w4.slider("약물화 가능성", 0, 100, 10),
        }
    ranked = rank_targets(rows, weights)
    st.session_state.ranked_targets = ranked
    target_table = []
    for rank, item in enumerate(ranked, start=1):
        target = item["target"]
        target_table.append(
            {
                "순위": rank,
                "유전자": target.get("approvedSymbol", ""),
                "단백질": target.get("approvedName", ""),
                "Ensembl ID": target.get("id", ""),
                "추천 점수 (0–1)": round(item["priority"], 4),
                "OT 연관": round(item["association"], 4),
                "유전학": round(item["genetics"], 4),
                "기존 약물": round(item["known_drug"], 4),
                "약물화 가능": "있음" if item["tractability"] else "자료 없음",
                "UniProt": item["uniprot"] or "연결 ID 없음",
            }
        )
    if target_table:
        st.dataframe(target_table, use_container_width=True, hide_index=True)
        st.caption("이 재순위 점수는 투명한 가중 합산 기준선이며, 학습된 머신러닝 모델이 아닙니다. Open Targets 상위 25개 안에서만 순위를 계산합니다.")
    else:
        st.info("이 질병에서 반환된 표적이 없습니다.")

    if ranked:
        labels = {
            f"{i + 1}. {item['target'].get('approvedSymbol', '?')} — {item['target'].get('approvedName', '')} ({item['uniprot'] or 'UniProt 연결 없음'})": item
            for i, item in enumerate(ranked)
        }
        chosen_target_label = st.selectbox("PDB 구조를 찾을 표적 단백질", list(labels))
        chosen_target = labels[chosen_target_label]
        st.markdown("### 표적 구조 찾기")
        if not chosen_target["uniprot"]:
            st.warning("선택한 표적에 연결된 UniProt accession이 Open Targets 응답에 없어 PDB 자동 검색을 할 수 없습니다.")
        elif st.button("이 단백질의 실험 PDB 구조 찾기"):
            try:
                with st.spinner("UniProt accession으로 RCSB PDB를 검색하고 구조 정보를 수집합니다..."):
                    uniprot_sequence = get_uniprot_sequence(chosen_target["uniprot"])
                    pdb_entities = search_pdb_by_uniprot(chosen_target["uniprot"])
                    pdb_structures = get_pdb_metadata(pdb_entities, chosen_target["uniprot"])
                    st.session_state.pdb_structures = attach_sequence_metrics(pdb_structures, uniprot_sequence)
                    st.session_state.pdb_target_id = chosen_target["target"].get("id")
            except Exception as exc:
                st.error(f"PDB 구조 검색 오류: {exc}")

    structures = st.session_state.get("pdb_structures", [])
    current_target_id = chosen_target["target"].get("id") if ranked else None
    if structures and st.session_state.get("pdb_target_id") == current_target_id:
        st.subheader("PDB 구조 비교")
        pdb_table = []
        for structure in structures:
            resolution = structure["best_resolution"]
            differences = structure["sequence_differences"]
            pdb_table.append(
                {
                    "PDB ID": structure["pdb_id"],
                    "구조 방법": structure["method_type"],
                    "해상도 (Å)": round(resolution, 2) if resolution is not None else "해당 없음",
                    "construct mutation 주석": structure["mutation_annotation"],
                    "UniProt 서열 커버리지": round(structure["reference_coverage"], 3) if structure["reference_coverage"] is not None else "계산 불가",
                    "기준서열 차이 후보": ", ".join(differences[:8]) if differences else "감지 안 됨/자료 없음",
                    "비단백질 성분 코드": ", ".join(structure["nonwater_components"][:12]) or "주요 성분 없음",
                    "구조 제목": structure["title"],
                }
            )
        st.dataframe(pdb_table, use_container_width=True, hide_index=True)
        st.caption("기준서열 차이는 UniProt canonical sequence와 PDB construct 서열의 국소 정렬 결과입니다. isoform·종 차이·construct 설계일 수도 있어 실제 mutation 확정 판정이 아닙니다. 해상도는 실험 방법 사이에서 단순 비교하지 마세요.")
        structure_labels = {}
        for structure in structures:
            resolution_label = f"{structure['best_resolution']:.2f} Å" if structure["best_resolution"] is not None else "해상도 미기재"
            structure_labels[f"{structure['pdb_id']} · {structure['method_type']} · {resolution_label}"] = structure
        selected_structure_label = st.selectbox("3D로 볼 구조 선택", list(structure_labels))
        selected_structure = structure_labels[selected_structure_label]
        st.markdown(f"#### {selected_structure['pdb_id']} — {selected_structure['title']}")
        st.caption("단백질은 cartoon, PDB에 포함된 비단백질 원자는 stick으로 표시합니다. 비단백질 성분이 모두 약물 리간드라는 뜻은 아닙니다.")
        try:
            render_structure(selected_structure["pdb_id"])
        except Exception as exc:
            st.error(f"3D 구조를 불러오지 못했습니다: {exc}")
        st.markdown(f"[RCSB PDB에서 구조 보기](https://www.rcsb.org/structure/{selected_structure['pdb_id']})")
    elif st.session_state.get("pdb_target_id") and ranked and st.session_state.get("pdb_target_id") != current_target_id:
        st.info("선택한 표적이 바뀌었어요. 해당 표적의 PDB 검색 버튼을 눌러 구조를 다시 불러와 주세요.")

    st.markdown("### MeSH 문헌 탐색")
    st.caption("선택 질병의 이름을 PubMed MeSH Terms 필드에서 검색합니다.")
    if st.button("질병 관련 논문 검색"):
        try:
            with st.spinner("PubMed를 검색 중입니다..."):
                st.session_state.pubmed_articles = search_pubmed_mesh(target_result.get("name", selected.get("name", "")))
        except Exception as exc:
            st.error(f"PubMed 검색 오류: {exc}")
    if "pubmed_articles" in st.session_state:
        articles = st.session_state.pubmed_articles
        if not articles:
            st.info("정확한 MeSH 검색 결과가 없어요. 질병의 동의어나 세부 MeSH 용어로 다시 검색할 수 있습니다.")
        for article in articles:
            authors = article["authors"] + (" et al." if article["authors"] else "")
            citation = f"{authors} {article['year']}. {article['journal']}. PMID: {article['pmid']}"
            if article["doi"]:
                citation += f" · DOI: {article['doi']}"
            st.markdown(f"- [{article['title']}]({article['url']})  \n  {citation}")

st.divider()
st.caption("KCD 표준용어 조회 · 질병/근거: Open Targets · 구조 검색/메타데이터: RCSB PDB 및 UniProt · 논문: NCBI PubMed")
