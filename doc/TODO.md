# TODO

- CG扩展到CH位点

>>Ignore Here
本页记录当前实现细节澄清、待决策事项和后续最小改造方向，优先服务 `mbias -> call` 阶段迭代。

## `methy_caller` 当前算法说明

- 当前 `scripts/methy_caller.py` 是一个 CpG-specific caller，不是通用 cytosine caller。
- 处理流程：
  1. 从参考基因组枚举指定染色体上的全部 `CpG` 位点。
  2. 按基因组坐标分 batch，对每个 batch 做 BAM pileup。
  3. 仅统计目标 flag `{99, 147, 83, 163}` 的 paired-end read。
  4. 跳过 deletion、refskip、缺失 query position、read 长度不足、落入 trimming 区间的观测。
  5. 对每个目标位点读取 `query_pos` 与 `query_pos + 1` 两个碱基，按二核苷酸分类：
     - `TG` -> methylated evidence
     - `CA` -> methylated evidence
     - `CG` -> unmethylated evidence
  6. 输出 headerless Bismark-like `.cov`：
     `chrom  pos  pos  methylation_percent  methylated_count  unmethylated_count`

## 为什么当前把 `TG/CA` 记作甲基化信号

- 当前脚本遵循 TAPS 语义，而不是 bisulfite 语义。
- 在 TAPS 中，甲基化 cytosine 最终读成 `T`，未甲基化 cytosine 保持 `C`。
- 对参考里的 `CpG` 位点：
  - 甲基化时，在一条链方向上可观察为 `TG`
  - 同一事件在另一条链/另一方向上可观察为 `CA`
  - 未甲基化时保持为 `CG`
- 因此当前实现采用：
  - methylated count = `TG + CA`
  - unmethylated count = `CG`

## `CH` calling 讨论结论

- 不能简单把当前 `CpG` 逻辑中的 `CG` 替换成 `CH`。
- 原因：`CpG` 具有双链对称性，允许把两条链折叠到 `TG/CA/CG` 这套固定二核苷酸规则里；`CH` 不具备这个性质。
- 如果要支持 `CH`，应改成 strand-aware 的 cytosine caller，而不是继续沿用 CpG-specific dinucleotide caller。

建议规则：

- 对正链 cytosine 位点（参考为 `C`）：
  - `C` = unmethylated
  - `T` = methylated
- 对反链 cytosine 位点（参考正向表现为 `G`，实质代表反链 C）：
  - `G` = unmethylated
  - `A` = methylated

## `CH` 支持的 MVP 路线

- 先做最小改造，不追求一步到位支持所有 context 折叠。
- 推荐 MVP：
  - 先只统计正链 `CH` 位点，即 `C[ACT]`
  - calling 时只看当前位置单碱基，不再依赖 `query_pos + 1` 的二核苷酸组合
  - 先输出 `strand=+`
  - 增加 context 标记，至少区分 `CA` / `CC` / `CT`，更好是 `CHG` / `CHH`
- 完整版再考虑：
  - 同时纳入反链 cytosine 位点
  - 按 cytosine event 统一汇总，而不是按固定二核苷酸汇总
  - 明确输出契约是否需要兼容 Bismark-like `.cov`，还是扩展为带 `strand/context` 的表

## 待办

- 确认 `call` 阶段后续是否仅继续支持 `CpG`，还是要引入 `CH`/`CHG`/`CHH`
- 若支持 `CH`，先确定 MVP 输出契约：是否新增 `strand` 与 `context` 列
- 评估是否将当前 caller 从 “CpG-specific dinucleotide caller” 重构为 “strand-aware cytosine caller”
- 若修改 calling 契约，同步更新 `README.md`、`TEST.md` 与下游 summary/文档说明
>>Ignore Here