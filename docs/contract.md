# 有限仿射位点数据

process(data)面向材料结构数据整理，接受显式operations和sites，展开[0,1)³分数坐标，输出可供后续结构处理使用的逐点来源和晶胞组成。不读取CIF或调用空间群库。

operations为1..24个唯一id，matrix为3×3整数（-2..2）且行列式±1，translation是三个有理数。每个有理数输入为[分子,分母]，分子±100000、分母1..100000，整数拒绝bool。平移与坐标按模1归一化；不允许同一变换用多个id重复。集合必须包含恒等、并在仿射复合模整数平移下闭合；不自动补全生成元。闭包、可逆矩阵与有限集共同确定本版运算群。

sites为1..200个唯一id，含species ASCII标记、point三个有理坐标、occupancy (0,1]。species是用户标签，不推断元素或化学价。相同源site的多个操作落到同一点，只计一次占位并保存所有operation id；不同源site即使同species，也保留独立贡献并相加。展开点总占位>1拒绝整个数据集，<=1允许混合占位和vacancy；用户不得把同一物理位点当作多个满占位不对称源重复输入。

输出按精确坐标元组排序，有理值为整数/n/d字符串。multiplicity为去重后轨道大小，stabilizer_size为操作数/轨道大小。cell_composition按全部展开贡献占位求和，不按操作次数重复计量。坐标严格相等才合并，不带容差。

这里只验证给定分数坐标仿射集合，不验证晶胞度量下是否保持距离，不自动判合法实际晶体、化学结构或衍射强度；没有动力学/粒子模拟。后续度量一致性、格式导入、衍射数据关联可接入已有operation/site/position来源。固定合成样例不代表真实材料。

输入严格字段校验；错误抛ValueError，无部分结果。io.publish只创建目标文件，已有目标和链接均不覆盖。固定example调用API，无通用CLI产品。

# 占位区间草稿与确值结构交付

solve_draft(data, delivery=None) 在 process 的操作/位点语义之上处理占位区间草稿：求出全晶胞约束下每条来源的精确单位范围；给出 delivery 分组请求时，判定能否交付一份可直接传给 process 的确值结构。全部离线、精确整数算术，不枚举可行分配，处理量不随刻度 K 的数值铺开。输入错误抛 ValueError；不修改调用方传入的数据。

## 草稿格式

- scale：整数刻度 K，1..10⁹（布尔值拒绝）。来源 s 分到 u_s 个单位时，其轨道每个点的占位都是 u_s/(mK)，m 为该轨道的互异展开点数。
- operations：与 process 相同的显式仿射操作、合法性、闭包与坐标语义。
- sites：1..200 条来源，字段 id/species/point 同 process；标量占位换成 units：[下界, 上界] 非空闭整数区间，端点 1..mK（mK 为该来源轨道的互异展开点数乘 K，展开与轨道划分沿用 process 语义：完整展开坐标集合相同即同轨道，同一点由多个操作产生只贡献一次）。所有来源都保留，不把缺测来源默认丢掉。
- orbit_observations：可选数组，每项 {"site": 成员来源ID, "units": [下界, 上界]}，约束该轨道全部来源单位之和；端点为 0..mK 的整数。引用同一轨道两次（即使经不同成员 ID）算重复，拒绝；未知来源拒绝。未列出的轨道仍受容量限制：来源单位之和最多为 mK。
- species_observations：可选数组，每项 {"species": 物种标签, "units": [下界, 上界]}，约束该物种的晶胞单位总数（该物种各来源单位之和）；端点为 0..K·M 的整数，M 为全草稿互异展开点数。同一物种重复声明、未知物种拒绝；未列出的物种无额外限制。
- 所有端点必须是整数，布尔值、浮点、字符串均拒绝；未知字段、空区间、越界端点均抛 ValueError。

## 求解结果

可行分配必须同时满足全草稿的所有区间与容量。响应 format 为 "occupancy-draft/1"：

- 无解（格式合法但约束矛盾）：{"format", "status": "infeasible"}，不输出范围或结构。
- 可行：status="feasible"，附 scale、cell_points（M）、orbits（按最小成员 ID 排序：index/multiplicity/sites）与 sources（按 ID 排序：id/species/orbit/units=[最小, 最大]）。每个来源的最小、最大单位数分别在某个完整可行分配中确实能达到，不是输入区间、局部估计或对单个示例的描述；两个端点不保证能彼此同时取到。

## 交付请求与结果

delivery 为 {"groups": [...]}（至少一组）；每组 {"members": [来源ID...], "representative": 组内来源ID}。交付可只选一部分轨道，但求解始终使用包括未选轨道在内的全部草稿与观测。对每个被选中的轨道，其全部来源必须被划分成互不重叠的非空组；同组只能收同物种且同轨道的来源。未知成员、重复归组、遗漏同轨道来源、跨物种或跨轨道合组、代表不在组内，均抛 ValueError，不默默调整划分。

- 对每组计算所有完整可行分配下组内单位和的实际最小、最大值（不按成员独立上下界分别相加）。交付的必要条件：每组单位和固定，且换成每点占位（固定值/(mK)）、约为最简分数后落在 process 有理数格式内（分子、分母 1..100000）；即使固定，超出旧格式范围也不舍入、不放宽旧接口。
- 满足条件时 delivery.status="delivered"：structure 为可直接传给 process 的完整输入（沿用原操作；每组只保留代表来源的 ID、坐标、物种，以固定的组占位替代区间）；groups 逐组给出 representative→members 映射、固定 units 与 occupancy，不让合并隐没来源。旧入口展开后，每个选中坐标上各物种的占位与任一原始可行分配在该坐标上的物种总占位一致，且不产生未选轨道的点。交付的是用户明确批准的投影，不是把未确定的原稿改写成唯一解。
- 任一组不固定或固定值超出旧格式时 delivery.status="blocked"，不给部分结构：blocking 逐组给出 representative/members/orbit/species、精确范围 units、原因（"not_fixed" 或 "occupancy_unrepresentable"，后者附约为最简分数的 occupancy）以及 witnesses——min/max 两个端点各自的一份完整可行分配（含未选中轨道的全部来源单位），可独立复算该组的端点总量。草稿可行而拒绝交付时，来源范围仍在 sources 中给出。

正常响应均可经 JSON 往返；同一输入仅重排操作、来源、观测、组或成员时，响应相同。

# 超胞分片与离线拼接

build_supercell_shard(data, basis, links, cells) 在 process 的原胞展开语义之上生成超胞分片包；assemble_supercell_shards(shards, finalize=False) 在接收端离线合并分片，不接受也不依赖原胞输入、原始路径、全胞清单或进程缓存。两处失败均抛ValueError，不修改调用方传入的数据，不返回半成品。

## 生成输入

- data：与 process 相同的 operations/sites 输入、校验与精确展开语义；特殊位点去重、混合占位求和、超占位拒绝在分片前生效，即使只请求一个胞片，超占位原胞也整体出错。
- basis：3×3 整数矩阵（JSON 数组），列向量为新基矢的旧晶胞分数坐标；每项绝对值≤10⁹，行列式 det 满足 1≤det≤10⁹。
- cells：1..64 个整数三元组，每项绝对值≤10⁹，为请求的旧原胞平移。两平移相差 basis 的整数组合即属同一胞片，同一请求内重复选中抛 ValueError。胞片对外一律使用规范代表 c：basis⁻¹c 各分量落在 [0,1) 的唯一整数向量；请求可用负数或非规范代表，输出只含规范代表。
- links：0..1000 条有向连接，字段 id/from/to/shift。id 遵循既有 ID 规则且唯一；from/to 为原胞输入的有理数坐标格式，模 1 后必须指向展开出的点；shift 为三维整数平移，每项绝对值≤10⁹。连接从旧胞 c 的 from 指向旧胞 c+shift 的 to；允许自环与同端点不同 id，不按端点去重，不做化学键推断。

## 包格式（format: "supercell-shard/1"）

- format：格式版本字符串；version：结构版本标识 "sha256:…"，由规范化的原胞展开数据（含逐点来源）、全部连接声明及基底决定，不受请求胞片、无语义的输入排序和等价坐标表示影响；不同版本不得混拼，即使边界记录仍完全吻合。
- basis、cell_count：基底与胞片总数（=det）。
- cells：本片覆盖的规范胞片代表，升序。
- points：每个选中胞片的全部展开点。像点身份为 (cell, point)，即规范胞片代表加原胞展开点坐标；position 为该点在新基底下各分量落入 [0,1) 的精确有理数；occupancy/vacancy 与 contributions（源位点 id、species、占位、产生它的操作 id）沿用 process 展开语义，不重复计数。
- links / out_links / in_links：完整连接实例与出端/入端记录。实例以 (link, cell)=（原连接 id、起点规范胞片代表）标识；记录含 from/to 两端身份（cell+point）与整向量 winding，满足 真实有向位移 = to 归一坐标 + winding − from 归一坐标；真实位移由原胞端点差及 shift 转到新基底得到，不从图上最短距离猜测，也不要求 winding 为零。两端都在选中胞片中只交付一个完整实例；只覆盖一端时分别交付可配对的出端或入端记录；跨片连接不删除，同胞片非零 winding 不算跨片。允许描述跨向未选胞片的连接，但不生成那些胞片的点清单。
- 数字格式：整数为 JSON 整数；有理数为 str(Fraction) 规范字符串（如 "3/4"、"0"）。points 按 (cell, point) 升序，连接记录按 (link, cell) 升序；同版本同覆盖下直接生成、乱序拼接、分阶段拼接产生等价且序列化确定的包。

## 拼接

- 输入为经 JSON 序列化再读回的包；逐项校验格式版本、基底行列式、规范代表、position 与 (cell, point) 一致性、占位守恒及本地引用：本地端点不存在、内边悬空、方向错误、两端均不属于本片等无效包均拒绝。
- 所有输入必须同 version 同 basis，覆盖互不重叠，合并最多覆盖 64 个胞片。
- 配对必须是同一连接实例的一出一入，两端及 winding 都吻合；另一端已在本次覆盖范围内却缺半边、内容不一致或重复均拒绝，不以位置相同代替配对。配完的只留下完整连接；仍越出覆盖范围的边界记录原样保留。无需验证两端记录被同时删除的情形，不要求签名或防伪。
- 普通拼接可只覆盖一部分，结果仍是同格式的包，可与后续分片继续拼接。
- finalize=True 只允许在覆盖基底行列式所代表的全部互异胞片且没有未配边界记录时返回完整周期结构，否则 ValueError；能拼出开放分片不等于能交付完整周期结构。
- 处理量与包大小不随未选中的胞片数或平行六面体的外接整数盒增长。

# 周期结构整组抽取与作业交付

extract_supercell_groups(package, repeats, requests, reserves=None, quota=None) 在 finalize 验收的超胞包之上，把包沿基底三个方向重复 nx×ny×nz 次，按请求点所属的连通组整组判定保留与额度、整组交付本次作业的仿真输入。全程离线、精确整数/有理数算术；不展开重复网格，不扫描无请求的连通组，不靠遍历超额组判断同组；只有确定交付的组才枚举其位点与连接（数量受额度约束），处理量不随重复网格体积或超额组体量增长。输入错误抛 ValueError，不修改调用方数据；正常的保留与额度不足是业务结果而非错误。

## 输入

- package：supercell-shard/1 超胞包，必须能通过 assemble_supercell_shards(..., finalize=True) 验收（覆盖齐全、无未配边界记录）；开放分片、缺覆盖包、篡改包、非法字段同样拒绝。包内记录顺序不计。
- repeats：[nx, ny, nz] 三个方向重复次数，各为 1..10⁹ 整数（布尔拒绝）。包内点（胞片代表 cell、原胞展开点坐标 point）加各分量落在 [0, n) 的重复像编号 replica 标识一个实际位点；replica 与 cell 不混用。
- requests：1..10000 条请求的非空数组，按首次出现顺序决定各组交付资格。每条 {"id", "cell", "point", "replica"}：id 遵循既有 ID 规则且列表内唯一；cell+point 必须引用包中已有的点（point 为包中规范有理数字符串）；replica 各分量为 [0, n) 整数。多条请求可指向同一实际位点。
- reserves：可选数组（缺省为空），0..10000 条保留标记，格式同请求；id 列表内唯一（与请求列表各自计数，可重名）；无顺序，重排不改变结果。
- quota：本次作业可接收的实际位点总数上限，整数 1..10000。

## 连通与归组

对包中每条完整有向连接，从任一重复像 r 出发有一条指向 (r+winding) mod n（逐分量取模）的边；连通组按无向可达定义，交付仍保留原方向、平行边与自环；连接只取包中声明，不按距离、物种或占位推断。同组的请求归并到结果的同一项（项按请求首次出现顺序排列），后续同组请求不再分配额度；即使该组最终不交付，也正确归组。

## 逐组判定（按请求首次出现顺序）

1. 组内任一实际位点带保留标记：status="reserved"，reserves 给出组内全部保留标记 ID（升序）；不输出位点或连接，不扣额度；即使组原本超额也以保留为准。标记不需要出现在请求列表里；同一包内点的某个重复像被标记不代表其所有重复像都被标记。
2. 否则仅当整组互异实际位点数不超过剩余额度时交付：status="delivered"，按实际位点数扣额度；不够则 status="quota_exceeded"，不输出部分组、不占额度，继续处理后面的组；未交付组不给出精确体量。

## 交付结果（format: supercell-extract/1）

顶层含 format、version（原结构版本）、basis（包基底）、repeats、groups。每个交付组自包含、可脱离原包读取：requests（组内请求 ID，按请求顺序）、status="delivered"、version、basis、repeats，以及：

- sites：组内全部实际位点，按 (cell, point, replica) 升序。每项含 cell/point/replica 身份、position=（包中 position + replica）逐分量除以重复次数的规范分数坐标（精确分数字符串），以及沿用原包的 occupancy/vacancy/contributions。
- links：组内全部连接实例，按 (link, cell, replica) 升序。实例身份为包内实例（连接 ID link、起点胞片代表 cell）加起点重复像 replica；from/to 为两个实际端点（cell+point+replica）；winding 为重复后基底中的整向量绕行量，满足 终点规范坐标 + winding − 起点规范坐标 = 包中原有向位移逐分量除以重复次数；不由最短距离猜测。取模后两端相同的边不合并。

保留组含 requests/status/reserves；超额组含 requests/status。结果可经 JSON 往返；包内记录与保留标记重排不改变结果，请求顺序具有业务含义。

# 周期迁移轨迹定位校准与时间窗交接

`handoff_windows(package, trajectory, windows)` 在 finalize 验收的超胞包（supercell-shard/1）之上解释一份压缩迁移轨迹，按调用方提供的时间窗独立核对两端实测像编号、判定测量冲突与区域越界，并把全部通过的窗口编译成接收方仅凭文档即可离线重建的有向路线。全程精确整数/有理数算术，不展开轨迹、不按重复次数逐轮推进、不扫描到首个越界步，也不为每个窗口重扫整张结点表。输入错误抛 `ValueError`，无部分结果、不修改调用方数据；测量冲突与越界是逐窗业务结果，互不影响。

## 轨迹格式

`trajectory = {"root": 根ID, "nodes": [结点...]}`。包内位点以 `(cell, point)` 标识，包基底下的规范分数坐标为该点记录的 `position`（记作 p）；轨迹状态另带不取模整数三元组像编号 r，实际分数坐标为 p+r。结点有三种：

- `{"id", "kind":"step", "link", "cell"}`：单跳，引用一条完整连接实例 `(link, cell)`。执行时从实例 `from` 位点到 `to` 位点，并把实例的 `winding` 加到 r。平行边不去重、不允许反向代用、不按最短距离重算绕行。
- `{"id", "kind":"sequence", "left", "right"}`：先左后右；左结点终点与右结点起点必须是同一包内位点身份（cell 与 point 均相同），否则拒绝整份轨迹。
- `{"id", "kind":"repeat", "child", "count"}`：连续执行 child 共 count 次，count 为 ≥2 的 JSON 整数（布尔拒绝）。要求 child 起止的包内位点身份相同，不要求净绕行为零。

结点 ID 遵循仓库 ID 规则且在结点表内唯一；定义无顺序、引用可向前；所有结点必须从根可达，不得有环（含自环）。上限：结点 ≤5000、根到叶路径结点数 ≤48、根实际步数 ≤10¹⁸，超限拒绝不截断；即使分支不在任何申请窗口内也完整校验。根不指定绝对起始像编号。

## 窗口输入

1..2000 个不同 ID 的窗口，每窗 `{"id", "l", "r", "left_image", "right_image", "region"}`：

- `l,r` 为根轨迹步号，JSON 整数，满足 0 ≤ l < r ≤ 总步数；步号 0 在第一跳前，步号 j 是走完第 j 跳后的状态。
- `left_image/right_image` 为该两步的实测像编号，三个 JSON 整数（布尔拒绝），每分量绝对值 ≤10⁶⁰，可为负。各窗测量独立，共享步号不要求相同像编号。
- `region` 为包基底下三轴闭区间 `[[lo,hi],...]`，端点为规范整数或有理数字符串（如 "0"、"-3/2"），不用浮点、无容差，lo ≤ hi。

## 判定顺序（逐窗独立）

1. **定位核对**：以左端实测像为该窗绝对起点，按轨迹压缩算术推得右端像。与右端实测不符即 `measurement_conflict`，回传 `measured_left/measured_right/derived_right`；不检查出界、不交接，绝不丢弃任一实测值后照常交接。
2. **区域判定**：两端相符后，对闭区间 [l,r] 内**所有步后状态**（含两端）做精确判定，实际坐标为该状态包内位点的 position 加绝对像编号。越界回 `out_of_region` 及**最早**越界状态的根轨迹步号 `step`、包内位点身份（cell/point/position）、实际像编号 `image` 与精确实际坐标 `coordinate`；不交付越界前缀。相同路线与净绕行换一对仍自洽的定位测量，可改变交接决定。
3. **交接**：全窗通过才为 `handoff`。路线从第 l 步校准状态出发，恰好执行第 l+1..r 跳，保留次序、原始连接实例与绕行，可推回右端状态。

## 交接文档（format: orbit-handoff/1）

顶层含 `format`、原结构 `version`、包 `basis`、`total_steps` 与按窗口 ID 排序的 `windows`（每项含 id/l/r 与状态及上述业务字段）。存在可交接窗口时另含 `handoff`，否则该键缺省（无路线数据）：

- `definitions`：跨窗共享、内容寻址（sha256 派生 ID）的压缩定义表，仅含获准路线实际需要的定义，按 ID 排序。三类：`hop`（引用 `(link, cell)` 实例）、`sequence`（left/right 两个定义 ID）、`repeat`（child 定义 ID 与整数 count，count ≥1）。不允许附上原轨迹未交接部分再让接收方跳过；窗口切片以共享子程序的部分定义表达，原结点子程序可跨窗复用而不按窗口份数复制。
- `links`：定义引用到的全部连接实例，含 `(link, cell)` 身份、两端身份与包内规范坐标 `position`、整向量 `winding`，按 `(link, cell)` 排序；同端点平行边各自保留。
- `routes`：按窗口 ID 排序，每条含原窗口步号 `l/r`、跳数 `hops`、根程序 ID `program`、起始位点身份与起始像编号 `image`、右端位点身份与右端像编号 `end_image`。

文档自包含：重建路线不需要原包、原轨迹或进程缓存。`rebuild_handoff(document)` 校验文档格式、引用闭包（所有定义从路线可达、无悬空或重复定义、hop 均有链接记录），返回的对象提供 `route_ids()`、`state_at(window_id, step)`（O(定义树深) 的状态查询）、`end_state(window_id)` 与惰性 `iter_hops(window_id)`（按序产出每跳的实例、winding 与两端绝对像编号，重复以计数展开、不预物化）。重建得到的末状态与文档 `end_image` 不一致时抛 `ValueError`。

正常结果可经 JSON 往返；仅重排包内记录、结点定义或窗口输入不改变输出。单窗定位与判定的工作量受图深与整数位长约束，不随总步数或重复轮数增长；交接数据规模同理。

# 本地语料发行与同源证据

`Corpus(path)` 打开（必要时创建）一个本地语料库文件，调用方持有路径，进程退出后可用同一路径重新打开。接口全部走 `ValueError` 子类 `CorpusError`；分折冲突另抛其子类 `FoldConflict`，异常属性 `conflicts` 给出冲突详情。布尔值不当作整数或 ID 接受。单写者语义：同一写入者进程内串行调用即可，不提供网络或多进程合并。

```python
from site_orbits import Corpus
with Corpus('corpus.db') as db:
    db.apply_batch({...})            # 样本、证据、撤回
    receipt = db.release('R', {'A': 'train', 'C': 'test'})
    manifest = db.fetch('R')          # 当前可取用清单
    audit = db.audit('R')             # 原清单 + 发行/撤下修订号
```

## 样本与几何等价

样本 ID 为既有 ID 字符规则内的唯一字符串；`data` 直接交给 `process(data)`，沿用其字段、展开、去重、占位合并和超占位/不闭合校验，不另设输入格式。每库最多 300 个样本，每样本最多 48 个展开位点；超限拒绝整批，不裁剪。

几何等价只看展开位点：把每点的来源贡献按物种相加得到精确占位向量，忽略源位点 ID、操作 ID、同物种占位拆写与输入顺序；存在同一个三维分数平移使两样本位点集合模 1 后逐点相同且对应点占位向量相同，即判等价。不忽略每点物种或占位，不用浮点容差，不把旋转、换轴、缩放、总组成相同或无连接当作等价，也不读晶格度量。

## 批次

`apply_batch(batch)` 接受对象，字段 `samples`、`evidence`、`retract` 均可省略，各自为数组：

- `samples: [{"id": 样本ID, "data": process输入}]`
- `evidence: [{"id": 证据ID, "a": 样本ID, "b": 样本ID}]`，两个端点必须不同；不同证据可指向同一对样本。
- `retract: [证据ID]`，只能撤回批前仍有效的证据；撤回只关闭该 ID，不删样本、不影响其他 ID、不取消几何等价。

成功批次返回 `{"revision": 新证据修订号, "samples": [...], "evidence": [...], "retracted": [...], "withdrawn": [...已撤下发行ID...]}`。修订号从空库的 0 起每成功批次加一。批次原子生效：空批、重复新 ID、样本 ID 已登记、证据 ID 已用过（一经分配不得复用）、悬空端点、撤回同批新证据或撤回无效证据、样本内容非法、任一上限超出，整批不生效；发行判断只用提交后的完整批次，不看中间态。同一批的新证据可引用同批新样本。

## 发行

`release(release_id, mapping)`：`release_id` 为唯一 ID；`mapping` 是非空映射，键为样本 ID、值为 `train`/`validation`/`test`，可用 `dict` 或 `[样本ID, 折]` 对列表，输入顺序不计。拒绝未知样本、重复样本、非法折、空映射。

判定：全库当前有效同源证据与几何等价合并为连通关系，连接路径可穿过本次未选中的样本；候选映射中同一连通组内被选中的样本必须同折。跨折则抛 `FoldConflict`，其 `conflicts[0]` 至少含 `source`、`target` 及完整 `steps` 证据链：每步 `{"from", "to", "kind"}`，几何步 `kind="geometry"` 并带精确分数平移 `translation`，同源步 `kind="evidence"` 并带证据 ID `evidence`。不自动移动折、不删样本、不以组号代替证据链。

成功发行返回回执并固定清单。回执字段：`format="corpus-receipt/1"`、`release_id`、`evidence_revision`（发行时修订号）、`sample_count`、`fold_counts`、`manifest_digest`（清单规范 JSON 的 sha256）。回执只证明曾经成功发行，不含可直接供训练读取的清单；清单由取用入口返回。

同一 `release_id` 配同一映射（忽略顺序）重试，一律先返回原回执，即使当前修订号或发行状态已变、原请求已不合规；不重查旧请求、不另写清单。同一 ID 配不同映射报错。失败发行不占用 ID。

## 取用、撤下与审计

- `fetch(release_id)` 返回该发行当前清单，格式 `corpus-manifest/1`：`release_id`、`evidence_revision`、`samples`（按样本 ID 升序，每项含 `sample`、`fold` 与展开位点 `positions`，来源贡献完整保留）。未知或已撤下的发行报错。
- 每次成功批次写入时，若某份仍可取用清单在新修订下首次出现跨折，该次提交在同一事务中把它撤下，并记录首次撤下修订号和**当时**的完整冲突证据链（含中间样本）；不影响未跨折的发行，即使其组扩张或合并。已撤下不能经 `fetch` 读取，后续证据撤回也不自动恢复；可用新发行 ID 提交现在合规分折。
- `audit(release_id)` 始终返回：`release_id`、`status`（`available`/`withdrawn`）、`receipt`、`manifest`（原清单）、`released_revision`、`withdrawn_revision`、`conflicts`（首次撤下依据或 null）。证据后来失效不改写撤下记录。不要求收回下游已取副本。
- 上限：每库最多 300 样本、2000 条累计证据 ID、100 次成功发行；超限拒绝。

## 持久化与重开

存储为单文件 SQLite（标准库），`meta`、`samples`、`evidence`、`releases` 四表。样本原样保存 `process` 展开结果与几何指纹；证据行、修订号、清单、回执和撤下记录在同一 `BEGIN IMMEDIATE` 事务内提交，提交前中断则回滚，未提交写入不可见。重新打开后不会出现“证据已提交但仍读到应撤下清单”或“有回执却无原清单”：证据修订号与撤下、回执与清单一起可见。提交后、调用方收到回执前中断的发行，按同 ID 重试规则返回已落下的一份回执与清单。

清单与回执均为排序确定的 JSON（`sort_keys=True`、紧凑分隔），可按字节序列化比较。
