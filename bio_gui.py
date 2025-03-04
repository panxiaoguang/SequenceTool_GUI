import dearpygui.dearpygui as dpg
from Bio.Seq import reverse_complement,translate
from Bio import SeqIO
import io
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import httpx
from Bio import Align
from Bio.Align import substitution_matrices

###################################################################################
#####################################################################################
def parse_coord(coord_string):
    chrom, position = coord_string.split(":")
    start, end = position.split("-")
    return chrom, int(start), int(end)

def get_genome_seq(genome, chrom, start, end, reverse) -> dict:
    base_url = "https://api.genome.ucsc.edu/getData/sequence"
    input_params = {
        "genome": genome,
        "chrom": chrom,
        "start": start,
        "end": end,
        "revComp": reverse
    }
    data = httpx.get(base_url, params=input_params).json()
    return data

def get_blat_data(userSeq, Type, db, output) -> dict:
    base_url = "https://genome.ucsc.edu/cgi-bin/hgBlat"
    input_params = {
        "userSeq": userSeq,
        "type": Type,
        "db": db,
        "output": output
    }
    data = httpx.get(base_url, params=input_params).json()
    return data

def global_alignment(ref, query, gap_open, gap_extend, end_gap_open, end_gap_extend, end_gap):
    matrix = substitution_matrices.load("NUC.4.4")
    if not end_gap:
        aligner = Align.PairwiseAligner(
            substitution_matrix=matrix,
            mode="global",
            internal_open_gap_score=0-gap_open,
            internal_extend_gap_score=0-gap_extend,
            end_gap_score=0
        )
    else:
        aligner = Align.PairwiseAligner(
            substitution_matrix=matrix,
            mode="global",
            internal_open_gap_score=0-gap_open,
            internal_extend_gap_score=0-gap_extend,
            end_open_gap_score=0-end_gap_open,
            end_extend_gap_score=0-end_gap_extend
        )
    ref = ref.strip().upper()
    query = query.strip().upper()
    alignments = aligner.align(ref, query)
    return alignments[0].format(), alignments[0].score

###################################################################################
#####################################################################################

# Create DearPyGui context
dpg.create_context()

# Create main window
dpg.create_viewport(title='SequenceTool', width=800, height=600)
dpg.set_viewport_pos([200, 100])
###############################################################
################功能回调区######################################
###############################################################
# Main window callback functions
def on_search_button_click():
    # 显示loading指示器
    dpg.configure_item("search_loading_indicator", show=True)
    sequence = dpg.get_value("search_sequence_input")
    genome_type = dpg.get_value("search_genome_selector")
    query_type = dpg.get_value("query_type_selector")
    
    newData = {}
    search_results = []

    genome_dict = {
        "Human": "hg38",
        "Zebrafish": "danRer11",
        "Mouse": "mm10"
    }
    newData['db'] = genome_dict[genome_type]
    newData['Type'] = query_type
    newData['userSeq'] = sequence
    newData['output'] = "json"
    
    blat_data = get_blat_data(**newData)
    raw_blat_data = blat_data['blat']
    if query_type == "protein":
        scale_factor = 3
    else:
        scale_factor = 1
    if raw_blat_data:
        for raw_data in raw_blat_data:
            matches = int(raw_data[0])
            repMatches = int(raw_data[2])
            misMatches = int(raw_data[1])
            qinsert = int(raw_data[4])
            tinsert = int(raw_data[6])
            score = scale_factor*(matches+round(repMatches/2)) - \
                scale_factor*misMatches-qinsert-tinsert
            search_results.append([str(score), str(int(raw_data[11])+1), raw_data[12], raw_data[8],
                                raw_data[13], str(int(raw_data[15])+1), raw_data[16], raw_data[18]])
    
    # 显示表格并添加新的搜索结果
    if search_results:
        # 隐藏文本显示区域，准备使用表格显示结果
        dpg.configure_item("search_result_text", show=False)
    
         # 清空现有表格内容
        for child in dpg.get_item_children("search_result_table", slot=1):
            dpg.delete_item(child)

        dpg.configure_item("search_result_table", show=True)
        for result in search_results:
            with dpg.table_row(parent="search_result_table"):
                for value in result:
                    dpg.add_text(value)
    else:
        dpg.configure_item("search_result_table", show=False)
    # 隐藏loading指示器
    dpg.configure_item("search_loading_indicator", show=False)

def on_fetch_button_click():
    # 显示loading指示器
    dpg.configure_item("loading_indicator", show=True)
    genome_type = dpg.get_value("genome_selector")
    coordinates = dpg.get_value("coordinates_input")
    direction = dpg.get_value("direction_radio")
    genome_dict = {
        "Human": "hg38",
        "Zebrafish": "danRer11",
        "Mouse": "mm10"
    }
    newData = {}
    newData['genome'] = genome_dict[genome_type]
    chrom, start, end = parse_coord(coordinates)
    newData['chrom'] = chrom
    newData['start'] = start-1
    newData['end'] = end
    if direction == "3'-5'":
        newData['reverse'] = 1
    else:
        newData['reverse'] = 0
    data = get_genome_seq(**newData)
    result = f">{data['chrom']}:{data['start']}-{data['end']}\n{data['dna']}"
    dpg.set_value("result_text", result)
    # 隐藏loading指示器
    dpg.configure_item("loading_indicator", show=False)

def on_reverse_complement_button_click():
    sequence = dpg.get_value("modify_sequence_input")
    out_list = []
    if ">" in sequence:
        for record in SeqIO.parse(io.StringIO(sequence), "fasta"):
            out_list.append(">"+str(record.id))
            out_list.append(reverse_complement(str(record.seq)))
        newseq = "\n".join(out_list)
    else:
        newseq = reverse_complement(sequence)
    dpg.set_value("modify_sequence_output", newseq)

def on_translate_button_click():
    sequence = dpg.get_value("modify_sequence_input")
    out_list = []
    if ">" in sequence:
        for record in SeqIO.parse(io.StringIO(sequence), "fasta"):
            out_list.append(">"+str(record.id))
            out_list.append(translate(str(record.seq)))
        newseq = "\n".join(out_list)
    else:
        newseq = str(translate(sequence))
    dpg.set_value("modify_sequence_output", newseq)

def on_upper_button_click():
    sequence = dpg.get_value("modify_sequence_input")
    out_list = []
    if ">" in sequence:
        for record in SeqIO.parse(io.StringIO(sequence), "fasta"):
            out_list.append(">"+str(record.id))
            out_list.append(str(record.seq).upper())
        newseq = "\n".join(out_list)
    else:
        newseq = sequence.upper()
    dpg.set_value("modify_sequence_output", newseq)

def on_lower_button_click():
    sequence = dpg.get_value("modify_sequence_input")
    out_list = []
    if ">" in sequence:
        for record in SeqIO.parse(io.StringIO(sequence), "fasta"):
            out_list.append(">"+str(record.id))
            out_list.append(str(record.seq).lower())
        newseq = "\n".join(out_list)
    else:
        newseq = sequence.lower()
    dpg.set_value("modify_sequence_output", newseq)

def on_clear_button_click():
    dpg.set_value("modify_sequence_input", "")
    dpg.set_value("modify_sequence_output", "")

def on_alignment_submit_click():
    Seq1 = dpg.get_value("seq1_input")
    Seq2 = dpg.get_value("seq2_input")
    gap_open = dpg.get_value("gap_open_input")
    gap_extend = dpg.get_value("gap_extend_input")
    end_gap = dpg.get_value("end_gap_input")
    end_gap_open = dpg.get_value("end_gap_open_input")
    end_gap_extend = dpg.get_value("end_gap_extend_input")
    # 执行全局比对
    alignments,score = global_alignment(Seq1, Seq2, gap_open, gap_extend, end_gap_open, end_gap_extend, end_gap)
    # 显示比对结果
    dpg.set_value("alignment_result_text", alignments)
###############################################################
################功能回调区######################################
###############################################################

with dpg.font_registry():
    # first argument ids the path to the .ttf or .otf file
    default_font = dpg.add_font("asserts/MiSans-Light.ttf", 18)
    second_font = dpg.add_font("asserts/cour-2.ttf", 15)

# Create main window content
with dpg.window(tag="main_window", pos=(0, 0), no_move=True, width=800, height=600):
    dpg.bind_font(default_font)
    # Create tabs
    with dpg.tab_bar():
        # Fetch Sequence tab
        with dpg.tab(label="Fetch Sequence"):
            with dpg.group(horizontal=False):
                # Genome type selector
                dpg.add_text("Select Species:")
                dpg.add_combo(
                    items=["Human", "Zebrafish", "Mouse"],
                    default_value="Human",
                    width=200,
                    tag="genome_selector"
                )
                
                # Genome coordinates input
                dpg.add_text("Coordinates:")
                dpg.add_input_text(
                    hint="For example: chr1:1000-2000",
                    width=400,
                    tag="coordinates_input"
                )
                
                # Sequence direction selection
                dpg.add_text("Select Direction:")
                with dpg.group(horizontal=True):
                    dpg.add_radio_button(
                        items=["5'-3'", "3'-5'"],
                        default_value="5'-3'",
                        tag="direction_radio"
                    )
                
                # Fetch sequence button
                dpg.add_button(
                    label="Fetch",
                    callback=on_fetch_button_click,
                    width=200,
                    height=30
                )
                
                # 添加loading指示器
                with dpg.group(horizontal=True):
                    dpg.add_loading_indicator(show=False, tag="loading_indicator")
                
                # Result display area
                dpg.add_text("Fetched Results:")
                rst = dpg.add_input_text(
                    multiline=True,
                    readonly=True,
                    width=750,
                    height=300,
                    tag="result_text"
                )
        
        # Other tabs (content to be implemented later)
        with dpg.tab(label="Search Seq from Genome"):
            with dpg.group(horizontal=False):
                # 序列输入框
                dpg.add_text("Input Sequence:")
                dpg.add_input_text(
                    multiline=True,
                    width=750,
                    height=100,
                    tag="search_sequence_input"
                )
                
                # 基因组类型选择器
                dpg.add_text("Select Genome Type:")
                dpg.add_combo(
                    items=["Human", "Zebrafish", "Mouse"],
                    default_value="Human",
                    width=200,
                    tag="search_genome_selector"
                )
                
                # 查询类型选择器
                dpg.add_text("Query Type:")
                dpg.add_combo(
                    items=["DNA", "protein", "translated RNA", "translated DNA"],
                    default_value="DNA",
                    width=200,
                    tag="query_type_selector"
                )
                
                # 搜索按钮
                dpg.add_button(
                    label="Search",
                    callback=on_search_button_click,
                    width=200,
                    height=30
                )
                
                # 添加loading指示器
                with dpg.group(horizontal=True):
                    dpg.add_loading_indicator(show=False, tag="search_loading_indicator")
                
                # 结果显示区域
                dpg.add_text("Search Results:")
                with dpg.table(header_row=True, resizable=True, policy=dpg.mvTable_SizingStretchSame,
                   borders_outerH=True, borders_innerV=True, borders_outerV=True, tag="search_result_table", show=False):
                    dpg.add_table_column(label="Score")
                    dpg.add_table_column(label="First")
                    dpg.add_table_column(label="Last")
                    dpg.add_table_column(label="Strand")
                    dpg.add_table_column(label="Chrom")
                    dpg.add_table_column(label="Start")
                    dpg.add_table_column(label="End")
                    dpg.add_table_column(label="Span")
                dpg.add_input_text(
                    multiline=True,
                    readonly=True,
                    width=750,
                    height=300,
                    tag="search_result_text"
                )
        
        with dpg.tab(label="Modify Sequence"):
            with dpg.group(horizontal=False):
                # 序列输入区域
                dpg.add_text("Input Sequence:")
                dpg.add_input_text(
                    multiline=True,
                    width=750,
                    height=100,
                    tag="modify_sequence_input"
                )
                
                # 功能按钮组
                with dpg.group(horizontal=True):
                    dpg.add_button(
                        label="Reverse Complement",
                        callback=on_reverse_complement_button_click,
                        width=150,
                        height=30
                    )
                    dpg.add_button(
                        label="Translate",
                        callback=on_translate_button_click,
                        width=150,
                        height=30
                    )
                    dpg.add_button(
                        label="To Upper",
                        callback=on_upper_button_click,
                        width=150,
                        height=30
                    )
                    dpg.add_button(
                        label="To Lower",
                        callback=on_lower_button_click,
                        width=150,
                        height=30
                    )
                    dpg.add_button(
                        label="Clear",
                        callback=on_clear_button_click,
                        width=150,
                        height=30
                    )
                
                # 结果输出区域
                dpg.add_text("Modified Results:")
                dpg.add_input_text(
                    multiline=True,
                    readonly=True,
                    width=750,
                    height=300,
                    tag="modify_sequence_output"
                )
        
        with dpg.tab(label="DNA Pairwise Alignment"):
            with dpg.group(horizontal=False):
                # 序列1输入区域
                dpg.add_text("Sequence 1:")
                dpg.add_input_text(
                    multiline=True,
                    width=750,
                    height=100,
                    tag="seq1_input"
                )
                
                # 序列2输入区域
                dpg.add_text("Sequence 2:")
                dpg.add_input_text(
                    multiline=True,
                    width=750,
                    height=100,
                    tag="seq2_input"
                )
                
                # 参数设置区域
                with dpg.group():
                    # Gap open参数
                    with dpg.group(horizontal=True):
                        dpg.add_text("Gap open:")
                        dpg.add_input_float(
                            default_value=10.0,
                            width=120,
                            tag="gap_open_input"
                        )
                    
                        # Gap extend参数
                        dpg.add_text("Gap extend:")
                        dpg.add_input_float(
                            default_value=0.5,
                            width=120,
                            tag="gap_extend_input"
                        )
                    with dpg.group(horizontal=True):
                        # End gap选项
                        dpg.add_text("End gap:")
                        dpg.add_checkbox(
                            default_value=True,
                            tag="end_gap_input"
                        )

                        # End gap open参数
                        dpg.add_text("End gap open:")
                        dpg.add_input_float(
                            default_value=10.0,
                            width=120,
                            tag="end_gap_open_input"
                        )

                        # End gap extend参数
                        dpg.add_text("End gap extend:")
                        dpg.add_input_float(
                            default_value=0.5,
                            width=120,
                            tag="end_gap_extend_input"
                        )
                
                # 提交按钮
                dpg.add_button(
                    label="Submit",
                    callback=on_alignment_submit_click,
                    width=200,
                    height=30
                )
                
                # 结果输出区域
                dpg.add_text("Alignment Results:")
                
                b2 = dpg.add_input_text(
                    multiline=True,
                    readonly=True,
                    width=750,
                    height=300,
                    tag="alignment_result_text"
                )
                dpg.bind_item_font(b2, second_font)

# 设置主窗口为主要窗口
dpg.set_primary_window("main_window", True)

# 设置DearPyGui
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()